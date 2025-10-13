import os
import re

def rename_tally_folders(directory="."):
    """
    Renames folders like:
        tallies_ARC_U001.50kgm3_Li07.5_294K
    into:
        tallies_ARC_294K_Li07.5_U001.50kgm3
    by moving the temperature and lithium enrichment forward.
    """
    for folder in os.listdir(directory):
        if not folder.startswith("volume_") or not os.path.isdir(os.path.join(directory, folder)):
            continue

        # Match pattern components
        match = re.match(r"(volume_)([A-Za-z]+)_(U[\d.]+kgm3)_(Li[\d.]+)_([0-9]+K)", folder)
        if match:
            prefix, breeder, uranium, lithium, temp = match.groups()
            new_name = f"{prefix}{breeder}_{temp}_{lithium}_{uranium}"
            old_path = os.path.join(directory, folder)
            new_path = os.path.join(directory, new_name)
            
            if old_path != new_path:
                os.rename(old_path, new_path)
                print(f"Renamed: {folder} → {new_name}")
        else:
            print(f"Skipped (no match): {folder}")


def delete_target_files(root_dir, target_name="fertile_n-gamma.csv"):
    """
    Recursively delete all files named 'fertile_n-gamma.csv' under root_dir.
    Prints each deletion and a summary.
    """
    deleted = 0
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if target_name in filenames:
            full_path = os.path.join(dirpath, target_name)
            try:
                os.remove(full_path)
                print(f"Deleted: {full_path}")
                deleted += 1
            except Exception as e:
                print(f"Could not delete {full_path}: {e}")
    print(f"\nDeleted {deleted} file(s) named '{target_name}'.")


if __name__ == "__main__":
    # rename_tally_folders(".")
    delete_target_files(".")
