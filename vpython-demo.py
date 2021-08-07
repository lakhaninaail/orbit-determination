import vpython
vpython.scene.autoscale = False
ball = vpython.sphere(radius=0.5)
ball.velocity = vpython.vector(25, 5, 0)

wallR = vpython.box(pos=vpython.vector(6, 0, 0),
                    size=vpython.vector(0.2, 12, 12),
                    color=vpython.color.red)
wallL = vpython.box(pos=vpython.vector(-6, 0, 0),
                    size=vpython.vector(0.2, 12, 12),
                    color=vpython.color.green)

varr = vpython.arrow(pos=ball.pos, axis=ball.velocity, color=vpython.color.yellow)

deltat = 0.005
while True:
    vpython.rate(50)
    ball.pos = ball.pos + ball.velocity * deltat

    if ball.pos.x > wallR.pos.x:
        ball.velocity.x = -ball.velocity.x
    if ball.pos.x < wallL.pos.x:
        ball.velocity.x = -ball.velocity.x

    varr.pos = ball.pos
    varr.axis = ball.velocity