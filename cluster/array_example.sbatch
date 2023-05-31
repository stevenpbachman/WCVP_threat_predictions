#!/bin/bash
#SBATCH --job-name=array_example
#SBATCH --array=1-5
#SBATCH --output=%A_%a.out

# Number to calculate square
numbers=(2 4 6 8 10)

# Get the current task ID
task_id=${SLURM_ARRAY_TASK_ID}

# Get the number for the current task
number=${numbers[$task_id - 1]}

# Calculate the square
square=$(($number * $number))

# Print the result
echo "Square of $number is $square"