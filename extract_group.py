#!/usr/bin/env python3
"""Extract group 8 from submission.csv"""

import os
import pandas as pd
GROUP_NUMBER = os.environ.get('GROUP_NUMBER', '8').zfill(3)
df = pd.read_csv('submission.csv')
group_number = df[df['id'].str.startswith(f'{GROUP_NUMBER}_')].copy()

print(f"Group {GROUP_NUMBER} has {len(group_number)} trees")

# Save to file
output_file = f'submission_{GROUP_NUMBER}.csv'
group_number.to_csv(output_file, index=False)
print(f"Extracted group {GROUP_NUMBER} to {output_file}")

