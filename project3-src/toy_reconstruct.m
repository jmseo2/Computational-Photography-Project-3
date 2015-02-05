function im_out = toy_reconstruct(im)

[imh, imw, nb] = size(im);
im2var = zeros(imh, imw);
im2var(1:imh*imw) = 1:imh*imw;

A = sparse([], [], []);
b = [];

e = 1;
for y = 1:imh
    for x = 1:imw-1
        A(e, im2var(y, x)) = 1;
        A(e, im2var(y, x + 1)) = -1;
        b(e) = im(y, x) - im(y, x + 1);
        e = e + 1;
    end
end

for y = 1:imh-1
    for x = 1:imw
        A(e, im2var(y, x)) = 1;
        A(e, im2var(y + 1, x)) = -1;
        b(e) = im(y, x) - im(y + 1, x);
        e = e + 1;
    end
end

A(e, im2var(1, 1)) = 1;
b(e) = im(1, 1);

b = b.';
v = A \ b;

im_out = zeros(imh, imw);
im_out(1:imh*imw) = v(1:imh*imw);




