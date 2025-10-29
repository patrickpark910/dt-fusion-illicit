
import os
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import matplotlib.pyplot as plt

try:
    import openmc
except Exception as e:
    raise SystemExit("OpenMC is required to run this script. Please run it in an environment with OpenMC installed.\n"
                     f"Import error: {e}")

def find_statepoint(start_dir: str = ".", prefer_latest: bool = True) -> Optional[str]:
    """Find a statepoint*.h5 file under start_dir. If prefer_latest, pick the most recent by mtime."""
    candidates: List[Path] = []
    for root, _, files in os.walk(start_dir):
        for f in files:
            if f.startswith("statepoint") and f.endswith(".h5"):
                candidates.append(Path(root) / f)
    if not candidates:
        return None
    if prefer_latest:
        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return str(candidates[0])
    return str(candidates[0])

def get_mesh_tally(sp: "openmc.StatePoint", tally_name: Optional[str] = None, score: Optional[str] = None) -> "openmc.Tally":
    """Return a mesh tally by name or the first tally that has a MeshFilter and (optionally) the requested score."""
    if tally_name is not None:
        try:
            t = sp.get_tally(name=tally_name)
            # verify it has a MeshFilter
            if any(isinstance(f, openmc.MeshFilter) for f in t.filters):
                if (score is None) or (score in t.scores):
                    return t
        except Exception:
            pass  # fall through to scanning all tallies

    # Scan all tallies
    for t in sp.tallies.values():
        has_mesh = any(isinstance(f, openmc.MeshFilter) for f in t.filters)
        has_score = (score is None) or (score in t.scores)
        if has_mesh and has_score:
            return t
    raise RuntimeError("No mesh tally found in the statepoint with the requested constraints.")

def mesh_info_from_filter(tally: "openmc.Tally") -> Tuple["openmc.RegularMesh", Tuple[int,int,int]]:
    """Extract the RegularMesh and its (nx, ny, nz) from the tally's filters."""
    for f in tally.filters:
        if isinstance(f, openmc.MeshFilter):
            mesh = f.mesh
            # Support RegularMesh (most common). If it's a newer Mesh type, try to extract shape.
            if hasattr(mesh, "dimension"):
                nx, ny, nz = mesh.dimension
            else:
                raise RuntimeError("Unsupported mesh type without 'dimension' attribute.")
            return mesh, (nx, ny, nz)
    raise RuntimeError("MeshFilter not found in tally.")

def get_score_index(tally: "openmc.Tally", score: Optional[str]) -> int:
    """Return the index of the requested score, defaulting to the first score if None."""
    if score is None:
        return 0
    try:
        return tally.scores.index(score)
    except ValueError:
        raise RuntimeError(f"Requested score '{score}' not present in tally scores: {tally.scores}")

def bin_centers_from_edges(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[-1:1])

def get_axis_edges(mesh: "openmc.RegularMesh") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y, z bin edges for the mesh."""
    # RegularMesh may provide lower_left, width and dimension;
    # edges = lower + np.arange(n+1)*width
    if not hasattr(mesh, "lower_left") or not hasattr(mesh, "width") or not hasattr(mesh, "dimension"):
        raise RuntimeError("Mesh missing required attributes (lower_left, width, dimension).")
    llx, lly, llz = mesh.lower_left
    wx, wy, wz = mesh.width
    nx, ny, nz = mesh.dimension
    x_edges = llx + np.arange(nx+1) * wx
    y_edges = lly + np.arange(ny+1) * wy
    z_edges = llz + np.arange(nz+1) * wz
    return x_edges, y_edges, z_edges

def plot_mesh_1d_profile(statepoint_path: Optional[str] = None,
                         tally_name: Optional[str] = None,
                         score: Optional[str] = "flux",
                         reduce_xy: bool = True,
                         save_path: Optional[str] = None):
    """
    Load a statepoint, find a mesh tally, and plot a 1D profile.
    If nx, ny, nz are all >1:
      - if reduce_xy=True, sum over x,y to plot axial (z) profile.
      - else, take a central slice at (ix, iy) = (nx//2, ny//2) and plot along z.
    If exactly one axis varies, plot along that axis.
    """
    # Resolve statepoint path
    sp_path = statepoint_path or find_statepoint(".")
    if sp_path is None:
        raise SystemExit("No statepoint*.h5 found under the current directory. "
                         "Please provide statepoint_path or run this where the statepoint exists.")

    print(f"Using statepoint: {sp_path}")
    with openmc.StatePoint(sp_path) as sp:
        tally = get_mesh_tally(sp, tally_name=tally_name, score=score)
        mesh, (nx, ny, nz) = mesh_info_from_filter(tally)
        sidx = get_score_index(tally, score)

        # Extract data arrays
        mean = tally.mean
        std  = getattr(tally, "std_dev", None)

        # If multiple scores or nuclides are present, squeeze down to selected score
        # mean shape is (n_filter_bins, n_scores, n_nuclides) typically; handle len>1 carefully
        arr = mean
        # reduce scores
        if arr.ndim >= 2 and arr.shape[1] > 1:
            arr = arr[:, sidx:sidx+1]
        elif arr.ndim >= 2 and arr.shape[1] == 1:
            # keep as is
            pass
        # reduce nuclides if present
        if arr.ndim == 3 and arr.shape[2] > 1:
            arr = arr[:, 0:1, 0]  # first nuclide
        elif arr.ndim == 3 and arr.shape[2] == 1:
            arr = arr[:, 0, 0]
        elif arr.ndim == 2:
            arr = arr[:, 0]
        elif arr.ndim == 1:
            pass
        else:
            arr = np.squeeze(arr)

        # reshape to mesh
        try:
            cube = arr.reshape((nx, ny, nz), order="C")
        except Exception:
            # try Fortran order if mapping differs
            cube = arr.reshape((nx, ny, nz), order="F")

        # Standard deviation handling
        err = None
        if std is not None:
            serr = std
            if serr.ndim >= 2 and serr.shape[1] > 1:
                serr = serr[:, sidx:sidx+1]
            elif serr.ndim >= 2 and serr.shape[1] == 1:
                pass
            if serr.ndim == 3 and serr.shape[2] > 1:
                serr = serr[:, 0:1, 0]
            elif serr.ndim == 3 and serr.shape[2] == 1:
                serr = serr[:, 0, 0]
            elif serr.ndim == 2:
                serr = serr[:, 0]
            elif serr.ndim == 1:
                pass
            else:
                serr = np.squeeze(serr)
            try:
                errcube = serr.reshape((nx, ny, nz), order="C")
            except Exception:
                errcube = serr.reshape((nx, ny, nz), order="F")
        else:
            errcube = None

        # Get axis edges and centers
        x_edges, y_edges, z_edges = get_axis_edges(mesh)
        x_centers = 0.5*(x_edges[:-1] + x_edges[1:])
        y_centers = 0.5*(y_edges[:-1] + y_edges[1:])
        z_centers = 0.5*(z_edges[:-1] + z_edges[1:])

        # Decide which axis to plot
        varying = [(nx>1), (ny>1), (nz>1)]
        if sum(varying) == 0:
            raise RuntimeError("Mesh has no varying dimension; cannot plot 1D profile.")
        if sum(varying) == 1:
            # Only one axis varies; plot along that axis directly
            if nx>1:
                prof = cube[:,0,0]
                errs = errcube[:,0,0] if errcube is not None else None
                x = x_centers; xlabel = "x (cm)"
            elif ny>1:
                prof = cube[0,:,0]
                errs = errcube[0,:,0] if errcube is not None else None
                x = y_centers; xlabel = "y (cm)"
            else:  # nz>1
                prof = cube[0,0,:]
                errs = errcube[0,0,:] if errcube is not None else None
                x = z_centers; xlabel = "z (cm)"
        else:
            # Multiple axes vary; prefer axial profile by reducing x,y
            if reduce_xy:
                prof = cube.sum(axis=(0,1))  # sum over x,y
                errs = errcube.sum(axis=(0,1)) if errcube is not None else None
                x = z_centers; xlabel = "z (cm) [x,y summed]"
            else:
                ix, iy = nx//2, ny//2
                prof = cube[ix, iy, :]
                errs = errcube[ix, iy, :] if errcube is not None else None
                x = z_centers; xlabel = f"z (cm) [slice at ix={ix}, iy={iy}]"

        # Plot
        plt.figure(figsize=(7,4.2))
        if errs is not None:
            plt.errorbar(x, prof, yerr=errs, fmt='o', markersize=3, linewidth=1)
        else:
            plt.plot(x, prof, '-o', markersize=3, linewidth=1)
        ttl = tally.name if tally.name else "mesh tally"
        scr = score if score is not None else tally.scores[0] if tally.scores else "score"
        plt.title(f"{ttl} — {scr}")
        plt.xlabel(xlabel)
        plt.ylabel("Mean tally value")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Saved plot to {save_path}")
        plt.show()

if __name__ == "__main__":
    # Default usage: search for latest statepoint and first mesh tally with 'flux' score
    plot_mesh_1d_profile()
