# color_utils.R

# Specific Palettes

palette_geno <- c("WT" = "#eecc67",
                  "Het" = "#98782d")

palette_sex <- c("M" = "#4478AB",
                 "F" = "#ED6677")

palette_learner <- c("Learner" = "#547B80",
                     "Non_Learner" = "#D1D3D4")

region_two_colors <- c("Hypothalamus" = "#332288",
                      "Other" = "#BBBBBB")


region_muted_colors <- c("Amygdala/Hypothalamus" = "#88CCEE", "Hypothalamus" = "#332288", "Midbrain" = "#DDCC77", 
                  "Striatum/Pallidum" = "#117733", "Thalamus" = "#CC6677", "Zona Incerta" = "#882255")


palette_class <- c("#4477AA", "#CCBB44", "#228833", "#AA3377")


palette_dissector <- c("Simona" = "#013220",
                       "Din" = "#ddf2d1")
# Large Palettes
palette_90 <- c("#FFE5B4", "#E6B8AF", "#F5DEB3", "#D8BFD8", "#B0E0E6", "#FFDAB9",
                "#E0FFFF", "#FFDEAD", "#FFE4E1", "#FFFACD", "#AFEEEE", "#FFF0F5",
                "#F0E68C", "#E6E6FA", "#F08080", "#FFF5EE", "#90EE90", "#FAEBD7",
                "#D3D3D3", "#FFB6C1", "#FFA07A", "#20B2AA", "#87CEFA", "#778899",
                "#B0C4DE", "#FFFFE0", "#32CD32", "#FAF0E6", "#800000", "#66CDAA",
                "#0000CD", "#BA55D3", "#9370DB", "#3CB371", "#7B68EE", "#00FA9A",
                "#48D1CC", "#C71585", "#191970", "#F5FFFA", "#FFE4B5", "#FFE4C4",
                "#FFEBCD", "#8A2BE2", "#A52A2A", "#DEB887", "#5F9EA0", "#7FFF00",
                "#D2691E", "#FF7F50", "#6495ED", "#FFF8DC", "#DC143C", "#00FFFF",
                "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400", "#BDB76B",
                "#8B008B", "#556B2F", "#FF8C00", "#9932CC", "#8B0000", "#E9967A",
                "#8FBC8F", "#483D8B", "#2F4F4F", "#00CED1", "#9400D3", "#FF1493",
                "#00BFFF", "#696969", "#1E90FF", "#B22222", "#FFFAF0", "#228B22",
                "#DCDCDC", "#F8F8FF", "#FFD700", "#DAA520", "#808080", "#ADFF2F",
                "#F0FFF0", "#FF69B4", "#CD5C5C", "#4B0082", "#98FB98", "#FFDAB9",
                "#EEE8AA", "#DDA0DD")

palette_100 <- c("#696969", "#A9A9A9", "#DCDCDC", "#2F4F4F", "#556B2F", "#6B8E23",
                 "#A0522D", "#A52A2A", "#2E8B57", "#800000", "#191970", "#006400",
                 "#708090", "#808000", "#483D8B", "#5F9EA0", "#008000", "#3CB371",
                 "#BC8F8F", "#663399", "#B8860B", "#BDB76B", "#008B8B", "#CD853F",
                 "#4682B4", "#D2691E", "#9ACD32", "#20B2AA", "#CD5C5C", "#00008B",
                 "#4B0082", "#32CD32", "#DAA520", "#7F007F", "#8FBC8F", "#B03060",
                 "#66CDAA", "#9932CC", "#FF0000", "#FF4500", "#00CED1", "#FF8C00",
                 "#FFA500", "#FFD700", "#6A5ACD", "#FFFF00", "#C71585", "#0000CD",
                 "#7CFC00", "#DEB887", "#40E0D0", "#00FF00", "#BA55D3", "#00FA9A",
                 "#8A2BE2", "#00FF7F", "#4169E1", "#E9967A", "#DC143C", "#00FFFF",
                 "#00BFFF", "#F4A460", "#9370DB", "#0000FF", "#F08080", "#ADFF2F",
                 "#FF6347", "#D8BFD8", "#B0C4DE", "#FF7F50", "#FF00FF", "#1E90FF",
                 "#DB7093", "#F0E68C", "#FA8072", "#EEE8AA", "#FFFF54", "#6495ED",
                 "#DDA0DD", "#ADD8E6", "#87CEEB", "#FF1493", "#7B68EE", "#F5DEB3",
                 "#AFEEEE", "#EE82EE", "#98FB98", "#7FFFD4", "#FF69B4", "#FFB6C1")


palette_43 <- c("#C69DAF", "#B984A7", "#B2B2B2", "#8CAC94", "#78A5A6", "#D7A69E", 
                "#A996BA", "#BEB7A9", "#FFF6CC", "#A1A964", "#94B5C8", "#C1A387", 
                "#C5C5C5", "#B5DFD9", "#C4C4C4", "#B17478", "#9AD1B5", "#A892AA",
                "#9E9F8D", "#A7DFF3", "#CCA277", "#C79DCD", "#7F7F7F", "#B5EAD7", 
                "#FFC0A9", "#AFAFAF", "#FFFFF0", "#859532", "#9CB4D6", "#FFA69E",
                "#A8C796", "#B8AF97", "#A4999F", "#8CD88C", "#BEBFEE", "#F0EDED", 
                "#F4ECDC", "#F0EAD6", "#F9E8DE", "#F2DCD1", "#E8DEE4", "#E7E2DA", 
                "#FDEAE1")

palette_70 <- c("#C69DAF", "#B984A7", "#B2B2B2", "#8CAC94", "#78A5A6", "#D7A69E", 
                "#A996BA", "#BEB7A9", "#FFF6CC", "#A1A964", "#94B5C8", "#C1A387", 
                "#C5C5C5", "#B5DFD9", "#C4C4C4", "#B17478", "#8C534E", "#657DAF", 
                "#F3CEB1", "#9AD1B5", "#A892AA", "#9E9F8D", "#A7DFF3", "#CCA277", 
                "#C79DCD", "#7F7F7F", "#B5EAD7", "#FFC0A9", "#AFAFAF", "#CC7431",
                "#FFFFF0", "#859532", "#9CB4D6", "#FFA69E", "#A8C796", "#B8AF97",
                "#A4999F", "#8CD88C", "#ED7AAB", "#FFDAC3", "#92B6D5", "#DCBDFB", 
                "#A5BDD5", "#A18B7E", "#D4C4E8", "#C9C0AB", "#FFCCB2", "#6F9BC6", 
                "#AFE8D6", "#ADB7DB", "#98DBB7", "#BDBDBD", "#EAE8E8", "#DCCCCC", 
                "#C69CFF", "#F7D6E4", "#FFDAC1", "#D3E1D6", "#B2C3B3", "#C8E6E0",
                "#8CB6D0", "#DFE6E6", "#DFD0C7", "#BFC5CC", "#90AFC5", "#C789D6",
                "#ACEADC", "#C9BCE6", "#FFD2BE", "#C5C575")

nxph4_palette <- c("#77AADD", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00",
                   "#EEDD88", "#EE8866", "#FFAABB")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_complexheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  draw(x)
  dev.off()
}