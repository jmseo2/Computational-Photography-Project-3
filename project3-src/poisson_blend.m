function im_out = poisson_blend(im_s, mask_s, im_background)

[imh, imw, nb] = size(im_s);
num_vars = sum(sum(mask_s));
im2var = zeros(imh, imw);

[ys, xs] = find(mask_s > 0);
im_out = im_background;

for i = 1:num_vars
    y = ys(i);
    x = xs(i);
    im2var(y, x) = i;
end

dy = [0, 0, -1, 1];
dx = [-1, 1, 0, 0];

for d = 1:nb
    A = sparse([], [], []);
    b = [];
    
    e = 1;
    for i = 1:num_vars
        y = ys(i);
        x = xs(i);
        
        if mask_s(y, x + 1) == 1
           A(e, im2var(y, x)) = 1;
           A(e, im2var(y, x + 1)) = -1;  
           b(e) = im_s(y, x, d) - im_s(y, x + 1, d);
           e = e + 1;
        end
        
        if mask_s(y + 1, x) == 1 
           A(e, im2var(y, x)) = 1;
           A(e, im2var(y + 1, x)) = -1; 
           b(e) = im_s(y, x, d) - im_s(y + 1, x, d);
           e = e + 1;
        end
        
        for j = 1:4
            ny = y + dy(j);
            nx = x + dx(j);
            if mask_s(ny, nx) == 0
                A(e, im2var(y, x)) = 1;
                b(e) = im_background(ny, nx, d) + im_s(y, x, d) - im_s(ny, nx, d);
                e = e + 1;
            end
        end
    end
    b = b.';
    v = lscov(A, b);
    
    for i = 1:num_vars
        im_out(ys(i), xs(i), d) = v(i);
    end
end
