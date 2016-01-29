for i = 3:10:173
    string = ['Glass_IPHys_a' threedig(i) '.X = (2.*Glass_IPHys_a' threedig(i-3) '_00.X + 3.*Glass_IPHys_a' threedig(mod(i+2,180)) '_00.X)./5'];
    eval(string)
    string = ['Glass_IPHys_a' threedig(i) '.Y = (2.*Glass_IPHys_a' threedig(i-3) '_00.Y + 3.*Glass_IPHys_a' threedig(mod(i+2,180)) '_00.Y)./5'];
    eval(string)
    string = ['Glass_IPHys_a' threedig(i) '.Field = (2.*Glass_IPHys_a' threedig(i-3) '_00.Field + 3.*Glass_IPHys_a' threedig(mod(i+2,180)) '_00.Field)./5'];
    eval(string)
    string = ['Glass_IPHys_a' threedig(i) '.Angle = (2.*Glass_IPHys_a' threedig(i-3) '_00.Angle + 3.*Glass_IPHys_a' threedig(mod(i+2,180)) '_00.Angle)./5'];
    eval(string)
end