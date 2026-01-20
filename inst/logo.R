library(hexSticker)
imgurl <-
    system.file("man/figures/GeneScout_logo3.png", package = "GeneScout")
sticker(
    imgurl,
    package = "",
    p_size = 14,
    s_x = 1,
    s_y = 1,
    s_width = .64,
    s_height = .8,
    p_color = "black",
    h_fill = "white",
    h_color = "white",
    filename = "man/figures/logo.png"
)
