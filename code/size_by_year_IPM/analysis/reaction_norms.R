##  A little script to visualize the implicit reaction norms imposed by
##  by our random year effects on growth and survival.


# Graphical parameters ----------------------------------------------------

set_graph_pars <- function(ptype = "panel4") {
  mgp <- c(2.5, 1, 0); mar <- c(4, 4, 1.5, 1); oma <- c(0, 0, 0, 0)
  switch(ptype,
         panel4 = par(mfrow=c(2,2), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.lab =1.2),
         panel2 = par(mfrow=c(1,2), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.axis=0.85),
         panel1 = par(mfrow=c(1,1), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.axis=0.85))
}

add_panel_label <- function(ltype="a",cex=1.2) {
  text <- paste(LETTERS[letters==ltype], ")", sep="")
  mtext(text=text, side=3, adj=0,cex=cex)
}


# Start PDF ---------------------------------------------------------------

pdf(file = "../../../Manuscripts/SizebyClimate/figures/reaction_norms.pdf",
    width = 8.5,
    height = 4)
set_graph_pars("panel2")

# Growth ------------------------------------------------------------------

grow_deriv <- function(x) {
  exp(x)
}

x <- seq(-2,2, length.out = 1000)
plot(x, grow_deriv(x), 
     type = "l", 
     xlab = expression(italic("E")), 
     ylab = expression(italic("g''(E)")))
text(0.2,4,"g is concave-up", cex = 0.8)
add_panel_label("a")

# Survival ----------------------------------------------------------------

logit_deriv <- function(x) {
  exp(x)*(1- exp(x))/(1+exp(x))^3
}
x <- seq(-8,8, length.out = 1000)
u <- exp(x)/(1+exp(x)); 
plot(u, logit_deriv(x), 
     type = "l", ylim=c(-0.15,0.15),
     xlab = expression(italic("Predicted survival S(E)")), 
     ylab = expression(italic("f ''(E)")))
text(0.25,0.12,"f is concave-up", cex = 0.8)
text(0.75,-0.12,"f is concave-down", cex = 0.8)
abline(v=0.5,col="blue",lty=2); 
add_panel_label("b")


# End PDF -----------------------------------------------------------------

dev.off()
