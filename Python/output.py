import openmc
import os, sys
import numpy as np
import pandas as pd

from Python.parameters import *
from Python.utilities import *


# Scores subtracted from 'scatter' to derive inelastic scattering.
# inelastic = scatter - elastic - (n,2n) - (n,3n) - (n,na) - (n,np)
SUBTRACTION_SCORES = ['elastic', '(n,2n)', '(n,3n)', '(n,na)', '(n,np)']


def _compute_inelastic(df):
    """
    Compute inelastic scattering by subtraction:
    inelastic = scatter - elastic - (n,2n) - (n,3n) - (n,na) - (n,np)
    Adds rows with score='inelastic'. Propagates uncertainties in quadrature.
    Original rows are preserved.
    """
    if df is None:
        return None

    if 'scatter' not in df['score'].values:
        return df

    group_cols = [c for c in df.columns if c not in ('score', 'mean', 'std. dev.')]

    base = df[df['score'] == 'scatter'][group_cols + ['mean', 'std. dev.']].set_index(group_cols)
    inel_mean = base['mean'].copy()
    inel_var = (base['std. dev.'] ** 2).copy()

    for score in SUBTRACTION_SCORES:
        sub = df[df['score'] == score]
        if sub.empty:
            continue
        sub = sub[group_cols + ['mean', 'std. dev.']].set_index(group_cols)
        inel_mean -= sub['mean'].reindex(inel_mean.index, fill_value=0.0)
        inel_var += (sub['std. dev.'].reindex(inel_mean.index, fill_value=0.0)) ** 2

    result = inel_mean.reset_index()
    result['std. dev.'] = np.sqrt(inel_var.values)
    result['score'] = 'inelastic'

    return pd.concat([df, result], ignore_index=True)


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

        # Compute inelastic = scatter - elastic - (n,2n) - (n,3n) - (n,na) - (n,np)
        raw['fertile']      = _compute_inelastic(raw['fertile'])
        raw['fertile_spec'] = _compute_inelastic(raw['fertile_spec'])
        raw['Pb']           = _compute_inelastic(raw['Pb'])
        raw['Pb_spec']      = _compute_inelastic(raw['Pb_spec'])

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

        leakage_spec = current_spec_df[(current_spec_df['cell'] == 99) & (current_spec_df['cellfrom'] != 99)].copy()

        leakage_Ebin_df = (leakage_spec.groupby(['energy low [eV]', 'energy high [eV]', 'energy mid [eV]'])
                                       .agg(mean=('mean', 'sum'), std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum())))
                                       .reset_index() )
        leakage_Ebin_df.columns = ['energy low [eV]', 'energy high [eV]',
                                   'energy mid [eV]', 'mean', 'std. dev.']
        return leakage_Ebin_df


    # Energy-binned neutron balance 

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

        U235_nufis_Ebin  = sum_over_cells(raw['fertile_spec'], 'U235',  'nu-fission')
        U238_nufis_Ebin  = sum_over_cells(raw['fertile_spec'], 'U238',  'nu-fission')
        Th232_nufis_Ebin = sum_over_cells(raw['fertile_spec'], 'Th232', 'nu-fission')

        U235_n2n_Ebin    = sum_over_cells(raw['fertile_spec'], 'U235',  '(n,2n)')
        U238_n2n_Ebin    = sum_over_cells(raw['fertile_spec'], 'U238',  '(n,2n)')
        Th232_n2n_Ebin   = sum_over_cells(raw['fertile_spec'], 'Th232', '(n,2n)')

        Li6_nt_Ebin    = sum_over_cells(raw['Li_spec'], 'Li6', '(n,Xt)')
        Li7_nt_Ebin    = sum_over_cells(raw['Li_spec'], 'Li7', '(n,Xt)')

        Be9_n2n_Ebin = sum_over_cells(raw['Be_spec'], 'Be9', '(n,2n)') if raw['Be_spec'] is not None else None
        Pb_n2n_Ebin  = sum_over_cells(raw['Pb_spec'], ['Pb204','Pb206','Pb207','Pb208'], '(n,2n)') if raw['Pb_spec'] is not None else None

        # Fertile inelastic scattering (MT 4 -> 'inelastic') and elastic scattering ('elastic')
        U235_inel_Ebin  = sum_over_cells(raw['fertile_spec'], 'U235',  'inelastic')
        U238_inel_Ebin  = sum_over_cells(raw['fertile_spec'], 'U238',  'inelastic')
        Th232_inel_Ebin = sum_over_cells(raw['fertile_spec'], 'Th232', 'inelastic')

        U235_elas_Ebin  = sum_over_cells(raw['fertile_spec'], 'U235',  'elastic')
        U238_elas_Ebin  = sum_over_cells(raw['fertile_spec'], 'U238',  'elastic')
        Th232_elas_Ebin = sum_over_cells(raw['fertile_spec'], 'Th232', 'elastic')
        Li6_elas_Ebin   = sum_over_cells(raw['Li_spec'], 'Li6', 'elastic')
        Li7_elas_Ebin   = sum_over_cells(raw['Li_spec'], 'Li7', 'elastic')
        Be9_elas_Ebin   = sum_over_cells(raw['Be_spec'], 'Be9', 'elastic') if raw['Be_spec'] is not None else None
        Pb_elas_Ebin    = sum_over_cells(raw['Pb_spec'], ['Pb204','Pb206','Pb207','Pb208'], 'elastic') if raw['Pb_spec'] is not None else None

        rxns_df = U235_ng_Ebin[['energy low [eV]', 'energy high [eV]', 'energy mid [eV]']].copy()

        channels = { 'U235(n,g)':     U235_ng_Ebin,
                     'U238(n,g)':     U238_ng_Ebin,
                     'Th232(n,g)':    Th232_ng_Ebin,
                     'U235(n,fis)':   U235_fis_Ebin,
                     'U238(n,fis)':   U238_fis_Ebin,
                     'Th232(n,fis)':  Th232_fis_Ebin,
                     'U235(n,nufis)':  U235_nufis_Ebin,
                     'U238(n,nufis)':  U238_nufis_Ebin,
                     'Th232(n,nufis)': Th232_nufis_Ebin,
                     'Li6(n,t)':     Li6_nt_Ebin,
                     'Li7(n,t)':     Li7_nt_Ebin,
                     'Be9(n,2n)':    Be9_n2n_Ebin,
                     'Pb(n,2n)':     Pb_n2n_Ebin,
                     'U235(n,2n)':   U235_n2n_Ebin,
                     'U238(n,2n)':   U238_n2n_Ebin,
                     'Th232(n,2n)':  Th232_n2n_Ebin,
                     'U235(n,inel)':  U235_inel_Ebin,
                     'U238(n,inel)':  U238_inel_Ebin,
                     'Th232(n,inel)': Th232_inel_Ebin,
                     'U235(n,el)':   U235_elas_Ebin,
                     'U238(n,el)':   U238_elas_Ebin,
                     'Th232(n,el)':  Th232_elas_Ebin,
                     'Li6(n,el)':    Li6_elas_Ebin,
                     'Li7(n,el)':    Li7_elas_Ebin,
                     'Be9(n,el)':    Be9_elas_Ebin,
                     'Pb(n,el)':     Pb_elas_Ebin,
                   } 

        for name, ch_df in channels.items():
            if ch_df is not None:
                merged = rxns_df.merge(
                    ch_df[['energy low [eV]', 'energy high [eV]', 'mean', 'std. dev.']],
                    on=['energy low [eV]', 'energy high [eV]'],
                    how='left')
                rxns_df[name]          = merged['mean'].fillna(0.0)
                rxns_df[f'{name}_sd'] = merged['std. dev.'].fillna(0.0)
            else:
                rxns_df[name]          = 0.0
                rxns_df[f'{name}_sd'] = 0.0

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

        # Multiplier (n,2n) and elastic — write both, one is zero
        if raw['Be'] is not None:
            Be = raw['Be']
            be9_n2n = Be[(Be['nuclide']=='Be9') & (Be['score']=='(n,2n)')][['cell','mean','std. dev.']]
            be9_list, be9_err = be9_n2n['mean'].tolist(), be9_n2n['std. dev.'].tolist()
            be9_el = Be[(Be['nuclide']=='Be9') & (Be['score']=='elastic')][['cell','mean','std. dev.']]
            be9_list_el, be9_err_el = be9_el['mean'].tolist() if len(be9_el) else zeros, be9_el['std. dev.'].tolist() if len(be9_el) else zeros
            pb_list,  pb_err  = zeros, zeros
            pb_list_el, pb_err_el = zeros, zeros
        else:
            Pb = raw['Pb']
            pb_n2n = sum_over_nuclides(Pb, '(n,2n)')
            pb_list,  pb_err  = pb_n2n['mean'].tolist(), pb_n2n['std. dev.'].tolist()
            pb_el = sum_over_nuclides(Pb, 'elastic')
            pb_list_el, pb_err_el = pb_el['mean'].tolist(), pb_el['std. dev.'].tolist()
            be9_list, be9_err = zeros, zeros
            be9_list_el, be9_err_el = zeros, zeros

        # Fertile capture — write both U238 and Th232, one is zero
        u238_ng  = fertile[(fertile['nuclide']=='U238')  & (fertile['score']=='(n,gamma)')][['cell','mean','std. dev.']]
        th232_ng = fertile[(fertile['nuclide']=='Th232') & (fertile['score']=='(n,gamma)')][['cell','mean','std. dev.']]

        # Fertile inelastic scattering (MT 4 -> 'inelastic')
        u235_inel  = fertile[(fertile['nuclide']=='U235')  & (fertile['score']=='inelastic')][['cell','mean','std. dev.']]
        u238_inel  = fertile[(fertile['nuclide']=='U238')  & (fertile['score']=='inelastic')][['cell','mean','std. dev.']]
        th232_inel = fertile[(fertile['nuclide']=='Th232') & (fertile['score']=='inelastic')][['cell','mean','std. dev.']]

        u235_inel_list  = u235_inel['mean'].tolist()  if len(u235_inel)  else zeros
        u235_inel_err   = u235_inel['std. dev.'].tolist()  if len(u235_inel)  else zeros
        u238_inel_list  = u238_inel['mean'].tolist()  if len(u238_inel)  else zeros
        u238_inel_err   = u238_inel['std. dev.'].tolist()  if len(u238_inel)  else zeros
        th232_inel_list = th232_inel['mean'].tolist() if len(th232_inel) else zeros
        th232_inel_err  = th232_inel['std. dev.'].tolist() if len(th232_inel) else zeros

        # Elastic scattering
        u235_elas  = fertile[(fertile['nuclide']=='U235')  & (fertile['score']=='elastic')][['cell','mean','std. dev.']]
        u238_elas  = fertile[(fertile['nuclide']=='U238')  & (fertile['score']=='elastic')][['cell','mean','std. dev.']]
        th232_elas = fertile[(fertile['nuclide']=='Th232') & (fertile['score']=='elastic')][['cell','mean','std. dev.']]
        li6_elas = Li[(Li['nuclide']=='Li6') & (Li['score']=='elastic')][['cell','mean','std. dev.']]
        li7_elas = Li[(Li['nuclide']=='Li7') & (Li['score']=='elastic')][['cell','mean','std. dev.']]

        u235_elas_list  = u235_elas['mean'].tolist()       if len(u235_elas)  else zeros
        u235_elas_err   = u235_elas['std. dev.'].tolist()  if len(u235_elas)  else zeros
        u238_elas_list  = u238_elas['mean'].tolist()       if len(u238_elas)  else zeros
        u238_elas_err   = u238_elas['std. dev.'].tolist()  if len(u238_elas)  else zeros
        th232_elas_list = th232_elas['mean'].tolist()      if len(th232_elas) else zeros
        th232_elas_err  = th232_elas['std. dev.'].tolist() if len(th232_elas) else zeros

        # Total fission summed over all nuclides
        all_fis   = sum_over_nuclides(fertile, 'fission')
        all_nufis = sum_over_nuclides(fertile, 'nu-fission')

        # Heating MW
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
            'U235(n,inel)':        u235_inel_list,
            'U235(n,inel)_stdev':  u235_inel_err,
            'U238(n,inel)':        u238_inel_list,
            'U238(n,inel)_stdev':  u238_inel_err,
            'Th232(n,inel)':       th232_inel_list,
            'Th232(n,inel)_stdev': th232_inel_err,
            'U235(n,el)':          u235_elas_list,
            'U235(n,el)_stdev':    u235_elas_err,
            'U238(n,el)':          u238_elas_list,
            'U238(n,el)_stdev':    u238_elas_err,
            'Th232(n,el)':         th232_elas_list,
            'Th232(n,el)_stdev':   th232_elas_err,
            'Li6(n,el)':           li6_elas['mean'].tolist()       if len(li6_elas) else zeros,
            'Li6(n,el)_stdev':     li6_elas['std. dev.'].tolist()  if len(li6_elas) else zeros,
            'Li7(n,el)':           li7_elas['mean'].tolist()       if len(li7_elas) else zeros,
            'Li7(n,el)_stdev':     li7_elas['std. dev.'].tolist()  if len(li7_elas) else zeros,
            'Be9(n,el)':           be9_list_el,
            'Be9(n,el)_stdev':     be9_err_el,
            'Pb(n,el)':            pb_list_el,
            'Pb(n,el)_stdev':      pb_err_el,
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
            'tot(n,fis)_stdev':    all_fis['std. dev.'].tolist(),
            'tot(n,nufis)':        all_nufis['mean'].tolist(),
            'tot(n,nufis)_stdev':  all_nufis['std. dev.'].tolist(),
            'heating [MW]':        heat_list,
            'heating_stdev':       heat_err_list,
            'fisq [MW]':           fisq_list,
            'fisq_stdev':          fisq_err_list,
        })

        totals = df.sum(numeric_only=True)
        totals['cell'] = 'total'
        df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
        return df


# Keep separate from OutputMixin class --ppark 2026-05-22

def collate_tallies(blanket, fertile_isotope, breeder_enrich, temp_k, vol_m3):
    """
    Collates all the tallies for given [blanket, fertile isotope, temperature]
    across multiple fertile kg/m³ into combined CSVs in ./Figures/Data/.

    Args:
        blanket (str): one of ['FLiBe', 'DCLL', 'HCPB']
        fertile_isotope (str): one of ['U238', 'Th232']
        breeder_enrich (float): enrichment of lithium-6 in breeder
        temp_k (float): temperature of the system
        vol_m3 (float): [m³] volume of breeder
    """

    rows_all = []
    rxnsE_list, fluxE_list, leakE_list = [], [], []

    tally_folders = [x for x in os.listdir("./OpenMC/")
                     if x.startswith(f"tallies_{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}")
                     and x.split("_")[-3].startswith(fertile_isotope)]

    dst = f"./Figures/Data/{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}_{fertile_isotope}"

    for folder in tally_folders:

        # Extract fertile loading from folder name
        part = folder.split("_")[-2]
        fertile = float(part.replace("kgm3", ""))
        mt = fertile * vol_m3 / 1e3

        # File paths (must match output.py exports)
        tally_summary = f"./OpenMC/{folder}/tallies_summary.csv"
        tally_leak    = f"./OpenMC/{folder}/tallies_leakage.csv"
        tally_rxnsE   = f"./OpenMC/{folder}/Ebins_rxns.csv"
        tally_fluxE   = f"./OpenMC/{folder}/Ebins_flux.csv"
        tally_leakE   = f"./OpenMC/{folder}/Ebins_leakage.csv"

        # Collate tallies_summary.csv
        try:
            df = pd.read_csv(tally_summary)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'tallies_summary.csv' not found in {folder}, skipping...")
            continue

        try:
            df_leak = pd.read_csv(tally_leak)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'tallies_leakage.csv' not found in {folder}, skipping...")
            continue

        tot  = df[df['cell'] == 'total']
        leak = df_leak[df_leak['void_cell'] == 'total']

        rows_all.append({
            'filename':        folder,
            'fertile_kg/m3':   fertile,
            'fertile_mt':      mt,
            'Li6(n,t)':        tot['Li6(n,t)'].values[0],
            'Li6(n,t)_sd':     tot['Li6(n,t)_stdev'].values[0],
            'Li7(n,Xt)':       tot['Li7(n,t)'].values[0],
            'Li7(n,Xt)_sd':    tot['Li7(n,t)_stdev'].values[0],
            'Be9(n,2n)':       tot['Be9(n,2n)'].values[0],
            'Be9(n,2n)_sd':    tot['Be9(n,2n)_stdev'].values[0],
            'Pb(n,2n)':        tot['Pb(n,2n)'].values[0],
            'Pb(n,2n)_sd':     tot['Pb(n,2n)_stdev'].values[0],
            'U238(n,g)':       tot['U238(n,g)'].values[0],
            'U238(n,g)_sd':    tot['U238(n,g)_stdev'].values[0],
            'Th232(n,g)':      tot['Th232(n,g)'].values[0],
            'Th232(n,g)_sd':   tot['Th232(n,g)_stdev'].values[0],
            'U235(n,inel)':    tot['U235(n,inel)'].values[0],
            'U235(n,inel)_sd': tot['U235(n,inel)_stdev'].values[0],
            'U238(n,inel)':    tot['U238(n,inel)'].values[0],
            'U238(n,inel)_sd': tot['U238(n,inel)_stdev'].values[0],
            'Th232(n,inel)':    tot['Th232(n,inel)'].values[0],
            'Th232(n,inel)_sd': tot['Th232(n,inel)_stdev'].values[0],
            'U235(n,el)':      tot['U235(n,el)'].values[0]       if 'U235(n,el)' in tot.columns else 0.0,
            'U235(n,el)_sd':   tot['U235(n,el)_stdev'].values[0] if 'U235(n,el)_stdev' in tot.columns else 0.0,
            'U238(n,el)':      tot['U238(n,el)'].values[0]       if 'U238(n,el)' in tot.columns else 0.0,
            'U238(n,el)_sd':   tot['U238(n,el)_stdev'].values[0] if 'U238(n,el)_stdev' in tot.columns else 0.0,
            'Th232(n,el)':     tot['Th232(n,el)'].values[0]      if 'Th232(n,el)' in tot.columns else 0.0,
            'Th232(n,el)_sd':  tot['Th232(n,el)_stdev'].values[0] if 'Th232(n,el)_stdev' in tot.columns else 0.0,
            'Li6(n,el)':       tot['Li6(n,el)'].values[0]        if 'Li6(n,el)' in tot.columns else 0.0,
            'Li6(n,el)_sd':    tot['Li6(n,el)_stdev'].values[0]  if 'Li6(n,el)_stdev' in tot.columns else 0.0,
            'Li7(n,el)':       tot['Li7(n,el)'].values[0]        if 'Li7(n,el)' in tot.columns else 0.0,
            'Li7(n,el)_sd':    tot['Li7(n,el)_stdev'].values[0]  if 'Li7(n,el)_stdev' in tot.columns else 0.0,
            'Be9(n,el)':       tot['Be9(n,el)'].values[0]        if 'Be9(n,el)' in tot.columns else 0.0,
            'Be9(n,el)_sd':    tot['Be9(n,el)_stdev'].values[0]  if 'Be9(n,el)_stdev' in tot.columns else 0.0,
            'Pb(n,el)':        tot['Pb(n,el)'].values[0]         if 'Pb(n,el)' in tot.columns else 0.0,
            'Pb(n,el)_sd':     tot['Pb(n,el)_stdev'].values[0]   if 'Pb(n,el)_stdev' in tot.columns else 0.0,
            'tot(n,fis)':      tot['tot(n,fis)'].values[0],
            'tot(n,fis)_sd':   tot['tot(n,fis)_stdev'].values[0],
            'tot(n,nufis)':    tot['tot(n,nufis)'].values[0],
            'tot(n,nufis)_sd': tot['tot(n,nufis)_stdev'].values[0],
            'tbr':             tot['tbr'].values[0],
            'tbr_sd':          tot['tbr_stdev'].values[0],
            'Pu239_kg/yr':     tot['Pu239_kg/yr'].values[0],
            'Pu239_kg/yr_sd':  tot['Pu239_kg/yr_stdev'].values[0],
            'U233_kg/yr':      tot['U233_kg/yr'].values[0],
            'U233_kg/yr_sd':   tot['U233_kg/yr_stdev'].values[0],
            'leak [n/src-n]':  leak['leakage'].values[0],
            'leak_sd':         leak['leakage_stdev'].values[0],
            'heat [MW]':       tot['heating [MW]'].values[0],
            'heat_sd':         tot['heating_stdev'].values[0],
            'fisq [MW]':       tot['fisq [MW]'].values[0],
            'fisq_sd':         tot['fisq_stdev'].values[0],
        })

        # Collate Ebins_rxns.csv (already summed over cells)
        try:
            df_rxnsE = pd.read_csv(tally_rxnsE)
            df_rxnsE['filename']   = folder
            df_rxnsE['fertile_mt'] = mt
            df_rxnsE['br_vol_m3']  = vol_m3
            rxnsE_list.append(df_rxnsE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_rxns.csv' not found in {folder}, skipping...")

        # Collate Ebins_flux.csv (still per-cell, needs groupby) 
        try:
            df_fluxE = pd.read_csv(tally_fluxE)
            cols = ['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'mean']
            fluxE = (df_fluxE[cols]
                     .groupby('energy mid [eV]', as_index=False)
                     .agg(**{'energy low [eV]':  ('energy low [eV]', 'first'),
                             'energy high [eV]': ('energy high [eV]', 'first'),
                             'mean':             ('mean', 'sum')}))
            fluxE['filename']      = folder
            fluxE['fertile_kg/m3'] = fertile
            fluxE['fertile_mt']    = mt
            fluxE['br_vol_m3']     = vol_m3
            fluxE_list.append(fluxE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_flux.csv' not found in {folder}, skipping...")

        # Collate Ebins_leakage.csv (already summed over cells)
        try:
            df_leakE = pd.read_csv(tally_leakE)
            df_leakE['filename']      = folder
            df_leakE['fertile_kg/m3'] = fertile
            df_leakE['fertile_mt']    = mt
            df_leakE['br_vol_m3']     = vol_m3
            leakE_list.append(df_leakE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_leakage.csv' not found in {folder}, skipping...")

    # Export collated files
    if rows_all:
        df_all = pd.DataFrame(rows_all).sort_values(by='fertile_kg/m3', ascending=True)
        df_all.to_csv(f"{dst}_summary.csv", index=False)

    if rxnsE_list:
        pd.concat(rxnsE_list, ignore_index=True).to_csv(f"{dst}_rxns.csv", index=False)

    if fluxE_list:
        df_fluxE_collated = pd.concat(fluxE_list, ignore_index=True)
        df_fluxE_collated = df_fluxE_collated[['filename', 'fertile_kg/m3', 'fertile_mt', 'br_vol_m3',
                                               'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]
        df_fluxE_collated.to_csv(f"{dst}_flux.csv", index=False)

    if leakE_list:
        pd.concat(leakE_list, ignore_index=True).to_csv(f"{dst}_leak.csv", index=False)

    print(f"{C.GREEN}Comment.{C.END} Collated tallies for {blanket} at {temp_k} K to: {dst}")