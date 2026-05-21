import openmc
import os, sys
import numpy as np
import pandas as pd

from Python.parameters import *
from Python.utilities import *


class OutputMixin:
    """
    Mixin class for tally extraction and collation.

    Provides extract_tallies() and its private helpers. Mixed into
    the Reactor class via multiple inheritance so that all methods
    have access to self.path, self.blanket_name, self.cells, etc.

    Usage in reactor.py:
        from Python.output import OutputMixin

        class Reactor(OutputMixin, ABC):
            ...
    """

    def extract_tallies(self):
        """Orchestrator: load statepoint, process tallies, export CSVs."""

        sp = find_sp(self.path, self.n_cycles)
        print(f"Reading tallies...")

        raw = self._load_tally_dfs(sp)

        leakage_df      = self._process_leakage(raw['current'])
        leakage_Ebin_df = self._process_leakage_spectrum(sp)
        rxns_df         = self._build_energy_binned_rxns(raw)
        summary_df      = self._build_summary_df(raw)

        print(summary_df.iloc[:, :8])

        # Export
        summary_df.to_csv(     f'./{self.path}/tallies_summary.csv', index=False)
        leakage_df.to_csv(     f'./{self.path}/tallies_leakage.csv',         index=False)
        rxns_df.to_csv(        f'./{self.path}/Ebins_rxns.csv',      index=False)
        leakage_Ebin_df.to_csv(f'./{self.path}/Ebins_leakage.csv',   index=False)

        flux_Ebin_df = raw['flux_spec'][raw['flux_spec']['cell'].between(30, 39, inclusive='both')]
        flux_Ebin_df.to_csv(   f'./{self.path}/Ebins_flux.csv',      index=False)


    # Load raw tally DataFrames from statepoint

    def _load_tally_dfs(self, sp):
        """
        Load all tally DataFrames from the statepoint and return them
        as a dict.  Be/Pb are loaded conditionally by blanket type;
        the absent multiplier is set to None.
        """
        raw = {
            'flux':         sp.get_tally(name='flux').get_pandas_dataframe(),
            'flux_spec':    sp.get_tally(name='flux spectrum').get_pandas_dataframe(),
            'Li':           sp.get_tally(name='Li rxn rates').get_pandas_dataframe(),
            'Li_spec':      sp.get_tally(name='Li rxn rates spectrum').get_pandas_dataframe(),
            'fertile':      sp.get_tally(name='Fertile rxn rates').get_pandas_dataframe(),
            'fertile_spec': sp.get_tally(name='Fertile rxn rates spectrum').get_pandas_dataframe(),
            'current':      sp.get_tally(name='current').get_pandas_dataframe(),
            'heating':      sp.get_tally(name='heating').get_pandas_dataframe(),
            'fisq':         sp.get_tally(name='fission-q').get_pandas_dataframe(),
            'Be':           None,
            'Be_spec':      None,
            'Pb':           None,
            'Pb_spec':      None,
        }

        if self.blanket_name in ['DCLL']:
            raw['Pb']      = sp.get_tally(name='Pb rxn rates').get_pandas_dataframe()
            raw['Pb_spec'] = sp.get_tally(name='Pb rxn rates spectrum').get_pandas_dataframe()
        else:
            raw['Be']      = sp.get_tally(name='Be rxn rates').get_pandas_dataframe()
            raw['Be_spec'] = sp.get_tally(name='Be rxn rates spectrum').get_pandas_dataframe()

        # Add energy-bin midpoints to all spectrum DataFrames
        for key in ['flux_spec', 'fertile_spec', 'Li_spec', 'Be_spec', 'Pb_spec']:
            if raw[key] is not None:
                raw[key]['energy mid [eV]'] = (raw[key]['energy low [eV]'] + raw[key]['energy high [eV]']) / 2

        return raw


    # Process leakage current into void cells

    def _process_leakage(self, current_df):
        """
        Compute total and per-void-cell leakage currents from blanket
        material cells into void cells.
        """
        void_ids     = [cell.id for cell in self.cells if cell.fill is None and cell.id != 10]
        material_ids = [23, 22, 42, 44]

        void_current_df = current_df[
            current_df['cell'].isin(void_ids) &
            current_df['cellfrom'].isin(material_ids)
        ].copy()

        total_leakage     = void_current_df['mean'].sum()
        total_leakage_err = np.sqrt((void_current_df['std. dev.']**2).sum())

        leakage_rows = []
        for vid in void_ids:
            mask = void_current_df['cell'] == vid
            leak_val  = void_current_df.loc[mask, 'mean'].sum()
            leak_err  = np.sqrt((void_current_df.loc[mask, 'std. dev.']**2).sum())
            neighbors = void_current_df.loc[mask & (void_current_df['mean'] > 0), 'cellfrom'].unique().tolist()
            leakage_rows.append({'void_cell': vid, 'cellfrom': str(neighbors),
                                 'leakage': leak_val, 'leakage_stdev': leak_err})

        leakage_rows.append({'void_cell': 'total', 'cellfrom': '',
                             'leakage': total_leakage, 'leakage_stdev': total_leakage_err})

        leakage_df = pd.DataFrame(leakage_rows)
        leakage_df.insert(0, 'fertile_kgm3', self.fertile_kgm3)
        return leakage_df


    # Process leakage current spectrum

    def _process_leakage_spectrum(self, sp):
        """
        Energy-binned leakage current into the outermost void (cell 99),
        summed over all contributing material cells.
        """
        current_spec_df = sp.get_tally(name='current spectrum').get_pandas_dataframe()
        current_spec_df['energy mid [eV]'] = (
            current_spec_df['energy low [eV]'] + current_spec_df['energy high [eV]']) / 2

        leakage_spec = current_spec_df[
            (current_spec_df['cell'] == 99) &
            (current_spec_df['cellfrom'] != 99)
        ].copy()

        leakage_Ebin_df = (
            leakage_spec
            .groupby(['energy low [eV]', 'energy high [eV]', 'energy mid [eV]'])
            .agg(mean=('mean', 'sum'),
                 std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum())))
            .reset_index()
        )
        leakage_Ebin_df.columns = ['energy low [eV]', 'energy high [eV]',
                                    'energy mid [eV]', 'mean', 'std. dev.']
        return leakage_Ebin_df


    # Energy-binned neutron balance ─────────────────────────────

    def _build_energy_binned_rxns(self, raw):
        """
        Assemble the energy-binned reaction rate DataFrame for all
        neutron balance channels.  Channels whose multiplier material
        is absent (Be in DCLL, Pb in FLiBe/HCPB) are filled with zeros.
        """
        
        U235_ng_Ebin   = sum_over_cells(raw['fertile_spec'], 'U235',  '(n,gamma)')
        U238_ng_Ebin   = sum_over_cells(raw['fertile_spec'], 'U238',  '(n,gamma)')
        Th232_ng_Ebin  = sum_over_cells(raw['fertile_spec'], 'Th232', '(n,gamma)')
        
        U235_fis_Ebin  = sum_over_cells(raw['fertile_spec'], 'U235',  'fission')
        U238_fis_Ebin  = sum_over_cells(raw['fertile_spec'], 'U238',  'fission')
        Th232_fis_Ebin = sum_over_cells(raw['fertile_spec'], 'Th232', 'fission')

        Li6_nt_Ebin    = sum_over_cells(raw['Li_spec'], 'Li6', '(n,Xt)')
        Li7_nt_Ebin    = sum_over_cells(raw['Li_spec'], 'Li7', '(n,Xt)')

        Be9_n2n_Ebin = sum_over_cells(raw['Be_spec'], 'Be9', '(n,2n)') if raw['Be_spec'] is not None else None
        Pb_n2n_Ebin  = sum_over_cells(raw['Pb_spec'], ['Pb204','Pb206','Pb207','Pb208'], '(n,2n)') if raw['Pb_spec'] is not None else None

        rxns_df = U235_ng_Ebin[['energy low [eV]', 'energy high [eV]', 'energy mid [eV]']].copy()

        channels = { 'U235_ng':   U235_ng_Ebin,
                     'U238_ng':   U238_ng_Ebin,
                     'Th232_ng':  Th232_ng_Ebin,
                     'U235_fis':  U235_fis_Ebin,
                     'U238_fis':  U238_fis_Ebin,
                     'Th232_fis': Th232_fis_Ebin,
                     'Li6_nt':    Li6_nt_Ebin,
                     'Li7_nt':    Li7_nt_Ebin,
                     'Be9_n2n':   Be9_n2n_Ebin,
                     'Pb_n2n':    Pb_n2n_Ebin,
                   } 

        for name, ch_df in channels.items():
            if ch_df is not None:
                merged = rxns_df.merge(
                    ch_df[['energy low [eV]', 'energy high [eV]', 'mean', 'std_dev']],
                    on=['energy low [eV]', 'energy high [eV]'],
                    how='left')
                rxns_df[name]            = merged['mean'].fillna(0.0)
                rxns_df[f'{name}_stdev'] = merged['std_dev'].fillna(0.0)
            else:
                rxns_df[name]            = 0.0
                rxns_df[f'{name}_stdev'] = 0.0

        rxns_df.insert(0, 'fertile_kg/m3', self.fertile_kgm3)
        return rxns_df


    # Per-cell scalar summary
    
    def _build_summary_df(self, raw):
        """
        Aggregate cell-wise scalar reaction rates, heating, and fissile
        production into a single summary DataFrame with a totals row.
        All columns use explicit isotope names (no generic 'mult' or 'fertile').
        """
        flux    = raw['flux']
        Li      = raw['Li']
        fertile = raw['fertile']
        heating = raw['heating']
        fisq    = raw['fisq']

        n_cells = flux['cell'].nunique()
        zeros   = [0.0] * n_cells

        # Lithium
        Li6_nt = Li[(Li['nuclide']=='Li6') & (Li['score']=='(n,Xt)')][['cell','mean','std. dev.']]
        Li7_nt = Li[(Li['nuclide']=='Li7') & (Li['score']=='(n,Xt)')][['cell','mean','std. dev.']]

        # Multiplier (n,2n) — write both, one is zero
        if raw['Be'] is not None:
            Be = raw['Be']
            be9_n2n = Be[(Be['nuclide']=='Be9') & (Be['score']=='(n,2n)')][['cell','mean','std. dev.']]
            be9_list,     be9_err = be9_n2n['mean'].tolist(), be9_n2n['std. dev.'].tolist()
            pb_list,      pb_err  = zeros, zeros
        else:
            Pb = raw['Pb']
            pb_n2n = (Pb[Pb['score']=='(n,2n)']
                      .groupby('cell')
                      .agg(mean=('mean', 'sum'),
                           **{'std. dev.': ('std. dev.', lambda x: np.sqrt((x**2).sum()))})
                      .reset_index())
            pb_list,      pb_err  = pb_n2n['mean'].tolist(), pb_n2n['std. dev.'].tolist()
            be9_list,     be9_err = zeros, zeros

        # Fertile capture — write both U238 and Th232, one is zero
        u238_ng  = fertile[(fertile['nuclide']=='U238')  & (fertile['score']=='(n,gamma)')][['cell','mean','std. dev.']]
        th232_ng = fertile[(fertile['nuclide']=='Th232') & (fertile['score']=='(n,gamma)')][['cell','mean','std. dev.']]

        # Total fission summed over all nuclides
        all_fis = (fertile[fertile['score']=='fission']
                   .groupby('cell')
                   .agg(mean=('mean', 'sum'),
                        std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum())))
                   .reset_index()
                   .sort_values('cell'))

        # Heating → MW
        heat_list     = (heating['mean']      * NPS_FUS * EV_TO_MJ).tolist()
        heat_err_list = (heating['std. dev.'] * NPS_FUS * EV_TO_MJ).tolist()
        fisq_list     = (fisq['mean']         * NPS_FUS * EV_TO_MJ).tolist()
        fisq_err_list = (fisq['std. dev.']    * NPS_FUS * EV_TO_MJ).tolist()

        # Fissile production — write both Pu239 and U233
        pu_scaling  = NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3
        u233_scaling = NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3

        u238_ng_list  = u238_ng['mean'].tolist()
        u238_ng_err   = u238_ng['std. dev.'].tolist()
        th232_ng_list = th232_ng['mean'].tolist()
        th232_ng_err  = th232_ng['std. dev.'].tolist()

        # Assemble
        cell_ids = [str(x) for x in flux['cell'].unique().tolist()]

        df = pd.DataFrame({
            'cell':                cell_ids,
            'flux':                flux['mean'].tolist(),
            'tbr':                 [a+b for a,b in zip(Li6_nt['mean'].tolist(), Li7_nt['mean'].tolist())],
            'tbr_stdev':           [a+b for a,b in zip(Li6_nt['std. dev.'].tolist(), Li7_nt['std. dev.'].tolist())],
            'U238(n,g)':           u238_ng_list,
            'U238(n,g)_stdev':     u238_ng_err,
            'Th232(n,g)':          th232_ng_list,
            'Th232(n,g)_stdev':    th232_ng_err,
            'Pu239_kg/yr':         [x * pu_scaling   for x in u238_ng_list],
            'Pu239_kg/yr_stdev':   [x * pu_scaling   for x in u238_ng_err],
            'U233_kg/yr':          [x * u233_scaling  for x in th232_ng_list],
            'U233_kg/yr_stdev':    [x * u233_scaling  for x in th232_ng_err],
            'Li6(n,t)':            Li6_nt['mean'].tolist(),
            'Li6(n,t)_stdev':      Li6_nt['std. dev.'].tolist(),
            'Li7(n,t)':            Li7_nt['mean'].tolist(),
            'Li7(n,t)_stdev':      Li7_nt['std. dev.'].tolist(),
            'Be9(n,2n)':           be9_list,
            'Be9(n,2n)_stdev':     be9_err,
            'Pb(n,2n)':            pb_list,
            'Pb(n,2n)_stdev':      pb_err,
            'tot(n,fis)':          all_fis['mean'].tolist(),
            'tot(n,fis)_stdev':    all_fis['std_dev'].tolist(),
            'heating [MW]':        heat_list,
            'heating_stdev':       heat_err_list,
            'fisq [MW]':           fisq_list,
            'fisq_stdev':          fisq_err_list,
        })

        totals = df.sum(numeric_only=True)
        totals['cell'] = 'total'
        df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
        return df