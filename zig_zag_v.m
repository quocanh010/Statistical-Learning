

%compute dct
function f = zig_zag_v(x, index_v)
    f = x.'
    f = f(:)
    f(:) = f(index_v)
end