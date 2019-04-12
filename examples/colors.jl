function hexcolor(rgb::UInt32)
    r = ((0xff0000 & rgb) >> 16) / 255
    g = ((0x00ff00 & rgb) >>  8) / 255
    b = ((0x0000ff & rgb)      ) / 255
    Plots.RGBA(r, g, b)
end

# Values taken from http://www.toutes-les-couleurs.com/code-couleur-rvb.php
troyes = hexcolor(0xfefdf0)
#troyes = Plots.RGBA(0.9961, 0.9922, 0.9412)
frambo = hexcolor(0xc72c48)
#frambo = Plots.RGBA(0.7804, 0.1725, 0.2824)
lichen = hexcolor(0x85c17e)
#lichen = Plots.RGBA(0.5216, 0.7569, 0.4941)
canard = hexcolor(0x048b9a)
#canard = Plots.RGBA(0.0157, 0.5451, 0.6039)
aurore = hexcolor(0xffcb60)
# colors taken from https://previews.123rf.com/images/capacitorphoto/capacitorphoto1410/capacitorphoto141000191/32438941-Graph-Icon-color-set-illustration--Stock-Vector.jpg
re = hexcolor(0xf59297)
gre = hexcolor(0xcbdf80)
blu = hexcolor(0x2394ce)
ora = hexcolor(0xfacb95)
yel = hexcolor(0xf2f08b)
