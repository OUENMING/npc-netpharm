
import sys
from pymol import cmd
rec, lig, out, color = sys.argv[1:5]
cmd.reinitialize(); cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("antialias", 2)
cmd.load(rec, "rec"); cmd.load(lig, "lig"); cmd.remove("solvent"); cmd.hide("everything")
cmd.show("cartoon", "rec"); cmd.color("grey80", "rec"); cmd.set("cartoon_transparency", 0.55)
cmd.select("pocket", "byres (rec within 4.0 of lig)"); cmd.show("sticks", "pocket"); cmd.color("lightblue", "pocket and elem C"); cmd.set("stick_radius", 0.09, "pocket")
cmd.show("sticks", "lig"); cmd.color(color, "lig and elem C"); cmd.set("stick_radius", 0.22, "lig"); cmd.util.cnc("all")
cmd.show("sticks", "organic"); cmd.set("depth_cue", 0); cmd.set("ray_shadow", 0)
cmd.zoom("lig", 6); cmd.set("orthoscopic", 1); cmd.png(out, 1200, 900, dpi=300, ray=1)
