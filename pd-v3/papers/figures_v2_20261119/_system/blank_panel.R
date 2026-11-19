# blank_panel.R - a pending panel: the block outlined with a thin line and nothing else (no text).
# Usage: BLANK_ID=2g Rscript render_panel.R blank_panel.R <block> <out_dir>
PANEL<-list(id=Sys.getenv('BLANK_ID','blank'),blocks=BLOCKS,what='pending content (outlined blank block)')
make_panel<-function(b)blank_grob(b$w,b$h)
