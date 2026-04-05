import os
import re

def cleanup_zone_identifiers(root_directory='.'):
    """
    Recursively finds and deletes all :Zone.Identifier files 
    starting from the root_directory.
    """
    deleted_count = 0
    
    for root, dirs, files in os.walk(root_directory):
        for file in files:
            if ":Zone.Identifier" in file:
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                    deleted_count += 1
                except OSError as e:
                    print(f"Error deleting {file_path}: {e}")
                    
    print(f"\nCleanup complete. Total files removed: {deleted_count}")


def rename_folders(directory="./"):

    # Define the target pattern and the replacement
    # We look for "U000" and replace it with "U238_000"
    search_str = "hcpb_"
    replace_str = "hcpb_Li60.0_"
    count = 0

    # List everything in the directory
    for foldername in os.listdir(directory):
        # Construct the full path
        old_path = os.path.join(directory, foldername)
        
        # Check if it is a directory and contains our search string
        if os.path.isdir(old_path) and search_str in foldername:
            new_name = foldername.replace(search_str, replace_str)
            new_path = os.path.join(directory, new_name)
            
            print(f"Renaming: {foldername} -> {new_name}")
            os.rename(old_path, new_path)
            count += 1
            
    print(f"\nFinished. {count} folders renamed.")



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
    print('1')
    # rename_folders()
    # cleanup_zone_identifiers()
    rename_folders()

    print('2')
    # delete_target_files(".")
