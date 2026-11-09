function cfg = brt_axis_config(axis_name, stack_file)
    cfg = FilterConfig.axisSpec(axis_name);
    grid = load(stack_file, 'grid_min', 'grid_max', 'grid_N');
    cfg.grid_min = grid.grid_min(:);
    cfg.grid_max = grid.grid_max(:);
    cfg.grid_num = grid.grid_N(:);
end
