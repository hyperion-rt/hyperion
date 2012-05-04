def find_last_iteration(file_handle):
    max_iteration = 0
    for group_name in file_handle:
        if "iteration" in group_name:
            iteration = int(group_name.split('_')[1])
            max_iteration = max(iteration, max_iteration)
    return max_iteration
