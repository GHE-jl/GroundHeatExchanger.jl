using Colors, CairoMakie

function update_fig_theme()
    fig_theme = Theme(
        #palette = (color = fig_color(), ),
        fontsize = 16,
        figure_padding = 7,
        Axis = (xlabelfont = "Times", ylabelfont = "Times",
            xticklabelfont = "Times", yticklabelfont = "Times",
            xgridvisible = false, ygridvisible = false,
            xtickalign = 1, ytickalign = 1,
            xminorticksvisible = true, yminorticksvisible = true,
            xminortickalign = 1, yminortickalign = 1),
        Legend = (framecolor = :transparent, backgroundcolor = :transparent,
            patchsize = (30, 10),
            labelfont = "Times", labelsize = 14)
    )
    set_theme!(fig_theme)
    return nothing
end

function fig_color()
    """
        fig_color()

    Outputs a palette of predefined colors for figures.
    """
    col = [
        RGB(30 / 255, 144 / 255, 205 / 255),      # 1 - Navy Blue
        RGB(255 / 255, 140 / 255, 0 / 255),       # 2 - Dark Orange
        RGB(50 / 255, 157 / 255, 13 / 255),       # 3 - Kinda Forest Green
        RGB(205 / 255, 175 / 255, 0 / 255),       # 4 - Gold yellow
        RGB(168 / 255, 40 / 255, 168 / 255),      # 5 - Purple
        RGB(205 / 255, 0 / 255, 0 / 255),         # 6 - Medium Red
        RGB(0 / 255, 0 / 255, 205 / 255),         # 7 - Medium Blue
        RGB(0 / 255, 205 / 255, 0 / 255),         # 8 - Medium green
        RGB(105 / 255, 105 / 255, 105 / 255),     # 9 - DimGrey1
        RGB(60 / 255, 60 / 255, 60 / 255)         # 10 - DimGrey2
    ]
    return col
end