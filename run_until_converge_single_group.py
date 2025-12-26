import subprocess
import shutil
import hashlib
import os
import sys

from datetime import datetime, timedelta
version = 'v21'


def file_hash(filepath):
    """Calculate MD5 hash of a file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()

def run_tree_packer():
    """Run tree_packer with live output and return final score."""
    print(f"Running tree_packer_{version}.exe...")
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

def move_file(iteration):
    """Move submission file to submission.csv."""
    src = f'submission_{version}.csv'
    dst = r'submission.csv'
    dst2 = f'submission_{version}_{iteration}.csv'
    if os.path.exists(src):
        shutil.copy(src, dst)
        # shutil.copy(src, dst2)
        return True
    return False

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
            print(f"tree_packer_{version}.exe failed!")
            break
        
        # Check if submission file exists
        submission_file = f'submission_{version}.csv'
        if not os.path.exists(submission_file):
            print(f"{submission_file} not created!")
            break
        
        # Check if score changed
        if previous_score is not None and final_score is not None:
            if abs(final_score - previous_score) < 1e-9:
                print("\nNo score change detected. Convergence achieved!")
                break
        
        # Check if files are different
        submission_file = f'submission_{version}.csv'
        new_hash = file_hash(submission_file)
        if initial_hash == new_hash:
            print("\nNo changes detected. Convergence achieved!")
            break
        
        # Update previous score
        if final_score is not None:
            previous_score = final_score
        
        # Move file
        move_file(iteration)
        print(f"Files differ. Moved {submission_file} -> submission.csv")
        # os.system('./bp.exe submission.csv submission.csv')
    
    print(f"\nCompleted after {iteration} iteration(s)")
    # from pymsgbox import alert
    # alert(f'Completed after {iteration} iteration(s)', f'Tree Packer {version}')
if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '32'
    group_number = os.environ.get('GROUP_NUMBER', '31')
    for n in range(100,10000):
        for r in range(10,1000):
            cmd = f'GROUP_NUMBER={group_number} ./single_group_optimizer  -n {n} -r {r} -o submission_v21.csv'
            main()
    import claude_notify
