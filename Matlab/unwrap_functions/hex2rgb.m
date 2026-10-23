function rgb = hex2rgb(h)
    rgb = sscanf(h(2:end), '%2x%2x%2x', [1 3])/255;
end
