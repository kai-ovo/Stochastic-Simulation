import os

def create_path(path):
    # Check if the directory does not exist, then create it
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory '{path}' created.")
    else:
        print(f"Directory '{path}' already exists.")