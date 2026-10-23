function saveFig(doSave, figDir, name)
    if doSave
        exportgraphics(gcf, fullfile(figDir, [name '.png']), 'Resolution', 150);
    end
end
