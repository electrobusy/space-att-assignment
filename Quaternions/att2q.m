function q = att2q(att)

q = [
    sin(att(1)/2).*cos(att(2)/2).*cos(att(3)/2)-cos(att(1)/2).*sin(att(2)/2).*sin(att(3)/2);
    cos(att(1)/2).*sin(att(2)/2).*cos(att(3)/2)+sin(att(1)/2).*cos(att(2)/2).*sin(att(3)/2);
    cos(att(1)/2).*cos(att(2)/2).*sin(att(3)/2)-sin(att(1)/2).*sin(att(2)/2).*cos(att(3)/2);
    cos(att(1)/2).*cos(att(2)/2).*cos(att(3)/2)+sin(att(1)/2).*sin(att(2)/2).*sin(att(3)/2);
    ];

end