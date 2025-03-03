import os

def delete_lines_before_restart(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    new_lines = []
    skip_count = 0

    for i in range(len(lines)):
        if skip_count > 0:
            skip_count -= 1

            continue

        if 'RESTART' in lines[i]:
            # Calculate the start index for the lines to be skipped

            start_index = max(0, i - 124)
            # Skip the lines before the 'RESTART' line

            skip_count = i - start_index

            continue

        new_lines.append(lines[i])

    with open(output_file, 'w') as file:
        file.writelines(new_lines)

# Usage

input_file = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/all/hematite-charges-1.hirshfeld'
output_file = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/all/test'
lines_delete = 128

delete_lines_before_restart(input_file, output_file)
