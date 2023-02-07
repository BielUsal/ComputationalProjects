using Makie 
epsilon=1
E=1
omega=1
odeSol(x,y) = Point(omega y, epsilon*E*cos(x)) # x'(t) = -x, y'(t) = 2y
scene = Scene(resolution =(400,400))
streamplot!(scene, odeSol, -2..2, -2..2, colormap = :plasma, 
    gridsize= (32,32), arrow_size = 0.07)
save("odeField.png", scene)
