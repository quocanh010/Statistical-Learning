

%compute dct
function index_v = zig_zag_v(x)
    x_dct = dct2(x);
    ind = reshape(1:numel(x_dct), size(x_dct)); %# indices of elements
    ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
    ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
    ind(ind==0) = [];                           %# keep non-zero indices
    v = x_dct(ind);
    [~, index_v] = max(abs(v(2:end)));
    index_v = index_v + 1;
end