function save_lightning_map(file_name,up,down)
    columnNames = {'Start_loc', 'Peak', 't12', 't13', 't23', 'CosAlphaOpt', 'CosBetaOpt', 'Azimuth', 'Elevation', 'Rcorr', 't123'};
    opts = detectImportOptions(file_name);
    opts.VariableNames = columnNames;
    table = readtable(file_name, opts);

    logicalIndex = abs(table.t123) < 1 & abs(table.Rcorr) > 0.3 & table.Start_loc > up & table.Start_loc < down;
    filteredTable1 = table(logicalIndex, :);
    Start_loc = filteredTable1.Start_loc;
    colorValues = (Start_loc - 3e8) / 2e8;
    set(0,'DefaultFigureVisible','off');
    figure;
    scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 1, colorValues, 'filled');
    title('Azimuth vs Elevation');
    xlabel('Azimuth');
    xlim([0, 360]);
    xticks(0:40:360);
    ylabel('Elevation');
    ylim([-40, 100]);
    yticks(-40:20:100);
    colormap('hsv');
    colorbar;
    caxis([0, 1.5]);
    grid on;
    [~, base_name, ~] = fileparts(file_name);
    png_file_name = [base_name '.png'];
    saveas(gcf,png_file_name);
end