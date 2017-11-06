import visit as v

v.LaunchNowin()

# open database
v.OpenDatabase("expanded_tags.vtk")

# add plot
v.AddPlot("Contour", "ww_n")

# adjust settings
att = v.ContourAttributes()
att.SetScaling(1)
att.minFlag = True
att.maxFlag = True
att.min = 5e-9
att.max = 0.5
v.SetPlotOptions(att)
print(att)

# draw the plot
v.DrawPlots()

#save the plot
s = v.SaveWindowAttributes()
s.format = s.PNG
s.fileName = "test"
s.width, s.height = 1024,768
s.screenCapture = 0
v.SetSaveWindowAttributes(s)
name = v.SaveWindow()
print(name)

# export


