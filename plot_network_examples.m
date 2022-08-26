% Plot some example fracture networks:
BASE_DIR = 'ensemble_1\\';

for ith_real = 1:10
    frac_data_full = load([BASE_DIR 'real_' num2str(ith_real) '.txt']);
    num_main_segm = size(frac_data_full, 1);
    frac_data = zeros(num_main_segm, 6);
    frac_data(:, [3:6]) = frac_data_full(:, [1:4]);

    figure()
    plot([frac_data(:, 3), frac_data(:, 5)]', [frac_data(:, 4), frac_data(:, 6)]', 'LineWidth', 2, 'color', 'black')
    hold on
    plot([frac_data(:, 3), frac_data(:, 5)]', [frac_data(:, 4), frac_data(:, 6)]', '.', 'LineWidth', 1, 'color', 'red')
    hold off
    axis equal
    grid on
    xlabel('x-coor', 'Interpreter', 'latex')
    ylabel('y-coor', 'Interpreter', 'latex')
    title('Fracture dataset', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14)
end