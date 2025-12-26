import subprocess
import shutil
import hashlib
import os
import sys

from datetime import datetime, timedelta

def file_hash(filepath):
    """Calculate MD5 hash of a file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()

def run_tree_packer():
    """Run tree_packer with live output and return final score."""
    print(f"Running single_group_optimizer.exe...")
    sys.stdout.flush()
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        shell=True,
        bufsize=1,
        universal_newlines=True
    )
    
    final_score = None
    # Print output line by line in real-time and parse score
    for line in process.stdout:
        print(line, end='', flush=True)
        # Parse "Final:   X.XXXXXX" line
        if line.startswith('Final:'):
            try:
                final_score = float(line.split()[1])
            except (ValueError, IndexError):
                pass
    
    process.wait()
    return process.returncode == 0, final_score

def main():
    iteration = 0
    previous_score = None
    start_time = datetime.now()
    while True:
        elapsed_time = datetime.now() - start_time
        print(f"Time elapsed: {elapsed_time}")
        if elapsed_time > timedelta(hours=11):
            print("Time limit exceeded!")
            break
        iteration += 1
        print(f"\n{'='*50}")
        print(f"Iteration {iteration}")
        print(f"{'='*50}\n")
        
        # Check hash of current submission.csv before processing
        initial_hash = file_hash('submission.csv')
        
        # Run tree_packer
        success, final_score = run_tree_packer()
        if not success:
            print(f"single_group_optimizer.exe failed!")
            break
        
        
        # Check if files are different
        new_hash = file_hash('submission.csv')
        if initial_hash == new_hash:
            print("\nNo changes detected. Convergence achieved!")
            break
    
    print(f"\nCompleted after {iteration} iteration(s)")
    # from pymsgbox import alert
    # alert(f'Completed after {iteration} iteration(s)', f'Tree Packer {version}')
if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '32'
    group_number = os.environ.get('GROUP_NUMBER', '31')
    cmd = f'./single_group_optimizer'
    main()
    import notify
