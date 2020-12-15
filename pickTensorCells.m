function cells = pickTensorCells(istart, ni, nx, ny)
    cells = (istart : (istart + ni - 1));
    cells = repmat(cells, ny, 1);
    cells = bsxfun(@plus, nx*(0 : (ny - 1))', cells);
    cells = reshape(cells', [], 1);
end