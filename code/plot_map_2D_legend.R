# --------------------------------------------------------------------
# code for complicated 2D legend


#' add_2d_legend
#'
#' adds a 2-dimensional legend to a current plot. 
#' assumes diverging colours for rows and changing transparency for columns
#'
#' @prop_y the proportion of the y dimension that is taken up with box
#' @prop_x the proportion of the x dimension that is taken up with box
#' @loc_y the proportion through the y dimension of the lower limit of the box
#' @loc_x the proportion through the x dimension of the left hand limit of the box
#' @row_colours character vector of colours for each row in the box
#' @row_label1 label for the higher colours (assuming diverging)
#' @row_label2 label for the lower colours (assuming diverging)
#' @row_arrows logical whether to add diverging arrows for the row labels
#' @col_trans numeric vector of transparency values for the columns of the box
#' @col_label character string of the name of the transparency dimension
#' @col_arrows logical whether to include the column arrows
#' @bg the background colour for the box
#' @label_cex the size of the legend labels 
#'
#' @return adds a 2D legend to the current plot
#'


add_2d_legend <- function(prop_y=0.15, prop_x=0.15, loc_y=0.15, loc_x=0.05,
							row_colours, row_label1, row_label2, row_arrows=TRUE,
							col_trans, col_label, col_arrows=TRUE,
							bg="grey80", label_cex=0.7, label_col="black"){

	# jiggery pokery to get the correct location of the legend in the plot area
	xdim <- 1/prop_x
	xstart <- - xdim*loc_x
	xend <- xdim + xstart

	ydim <- 1/prop_y
	ystart <- - ydim*loc_y
	yend <- ydim + ystart

	# getting the divisions between cells
	nr <- length(row_colours)
	box_y_start_end <- (1/nr)*0:nr

	nc <- length(col_trans)
	box_x_start_end <- (1/nc)*0:nc

	# adding the grid of colours
	par(new=TRUE)
	plot(0, 0, col=alpha("white", 0), xlim=c(xstart, xend), ylim=c(ystart, yend), 
		xaxt="n", yaxt="n", bty="n", xlab="", ylab="")

	polygon(x=c(0, 1, 1, 0, 0), y=c(0, 0, 1, 1, 0), col=bg, border=alpha("white", 0))

	for(i in 1:length(row_colours)){
		ys <- c(box_y_start_end[i], box_y_start_end[i+1]) 

		for(j in 1:length(col_trans)){
			xs <- c(box_x_start_end[j], box_x_start_end[j+1])

			polygon(x=c(xs[1], xs[2], xs[2], xs[1], xs[1]), y=c(ys[1], ys[1], ys[2], ys[2], ys[1]),
				col=alpha(row_colours[i], alpha=col_trans[j]), border=alpha("white", 0))
		}
	}

	# these could be optimised to be proportional to the extent (so always a
	# given distance from the coloured legend box)

	# --------------------------------------------------------------------
	# adding the labels and arrows

	# figure out pretty places to put them
	arrow_loc_x <- 1 + 0.02*xdim
	text_loc_x <- 1 + 0.02*xdim

	arrow_loc_y <- 1 + 0.02*ydim
	text_loc_y <- 1 + 0.04*ydim

	# along the top
	if(col_arrows) arrows(x0=0, y0=arrow_loc_y, x1 = 1, length=0.05, col=label_col)
	text(x=0.5, y=text_loc_y, labels=col_label, cex=label_cex, col=label_col)

	# on the side (diverging arrows from the middle)
	if(row_arrows) arrows(x0=arrow_loc_x, y0=0.55, y1=1, length=0.05, col=label_col)
	text(x=text_loc_x, y=0.75, pos=4, labels=row_label1, cex=label_cex, col=label_col)

	if(row_arrows) arrows(x0=arrow_loc_x, y0=0.45, y1=0, length=0.05, col=label_col)
	text(x=text_loc_x, y=0.25, pos=4, labels=row_label2, cex=label_cex, col=label_col)

}




# ####################################################################
# EXAMPLE PLOT

# loop through each of the 10 transparency categories for(i in 1:10)

# "transparency_raster" is a raster where the values are the degree of 
# transparency desired in each location (0-1)

# "colour_hue_raster" is a raster where the values are the colour desired 
# in each location. So will be a scale from 1 - length (colour)
# Actually, I think the code could probably cope with any scale.  



# --------------------------------------------------------------------
# colour settings

col_scale <- rev(viridis(10))

xlimits <- c(-7000000, 6000000)
ylimits <- c(-6800000, 8800000)

xrange <- xlimits[2] - xlimits[1]
yrange <- ylimits[2] - ylimits[1]

w_cm <- 14
h_cm <- w_cm*(yrange/xrange)

div_col <- cm.colors(9)

col_land <- "grey70" # light grey land background
col_sea <- "grey20"  # dark grey sea
col_text <- "white"  # white text (needs to contrast with sea colour)

library(scales)



# --------------------------------------------------------------------
# plot pd with pi as a transparency modifier

pd_plot_loc <- paste0(output_folder, "Map_PDPI_", var, "_week", week_code, ".png")
png(pd_plot_loc, width=w_cm, height=h_cm, units="cm", pointsize=9, res=600)

	par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))

	# plot background points of N and S America
	image(template_plot_new, maxpixels=ncell(template_plotproj),
		col=c(col_land, col_sea), xlim=xlimits, ylim=ylimits, zlim=c(0,1),
		xaxt="n", yaxt="n", 
		axes=FALSE, box=FALSE, bg=col_sea)

	# loop through different values of PI and plot each with different transparency

	for(i in 1:10){
		min_pi <- i*0.1-0.1
		max_pi <- i*0.1

		this_band <- transparency_raster
		this_band[is.na(transparency_raster)] <- 0
		this_band[!is.na(transparency_raster)] <- 0
		this_band[transparency_raster > min_pi & transparency_raster <= max_pi] <- 1

		colour_transparency_raster <- colour_hue_raster * this_band
		colour_transparency_raster[colour_transparency_raster==0] <- NA

		div_col_trans <- alpha(div_col, alpha=max_pi)

		image(colour_transparency_raster, maxpixels=ncell(template_plotproj),
			col=div_col_trans, xlim=xlimits, ylim=ylimits, zlim=c(0, 1), 
			xaxt="n", yaxt="n", 
			axes=FALSE, box=FALSE, # alpha=max_pi,
			add=TRUE) 

	} # close i


	# add title 
	text(x=xlimits[1], 
		y=ylimits[1]+ 0.12*(ylimits[2]-ylimits[1]),
		labels=paste("PD x PI", var), pos=4, font=1, cex=0.5, col=col_text)

	text(x=xlimits[1], 
		y=ylimits[1]+ 0.09*(ylimits[2]-ylimits[1]),
		labels=paste("week", week), pos=4, font=1, cex=0.5, col=col_text)

	# add complicated 2D legend. 

	add_2d_legend(prop_y=0.15, prop_x=0.15, loc_y=0.20, loc_x=0.05,
		row_colours=div_col, row_label1="positive effect", row_label2="negative effect", row_arrows=TRUE,
		col_trans=alpha.val, col_label="predictor importance", col_arrows=TRUE,
		bg=col_land, label_cex=0.6, label_col=col_text)


dev.off()






