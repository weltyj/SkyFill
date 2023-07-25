#Lens Flare
Lensflare is a small utility added as an option to skyfill generate lens flare artifacts caused by the sun or other bright lights refracting/reflecting through camera lenses.

It is not an attempt to create a true physical model of how flares are created in a lens, it is just a simple way to generate some flare and aperture ghosts  that is good enough to have a pleasantly realistic rendering

Lensflare reads data from an external file, and applies the results to a current image buffer that is open in skyfill.

##Datafile Format
The data file is a simple text file, with a small set of "commands",
The commands are briefly

- FILE - insert lines from an additional file
- BLACKSKY -set background image to black for testing
- LEAVES -set number of leaves on camera aperture
- FLARE -set ray types, from simple to more complicated
- SCALE -scale the dimensions of rays or ghosts that follow in the data file)
- ROTATE -rotate about the image center the rays or ghosts that follow in the data file)
- IMAGE_CENTER - set the location of the image center to form the flares/ghosts
- GHOST_LINE - create a set of aperture ghosts in a line extending from the image center
- GHOST_RING - create a set of aperture ghoss in a ring around the image center

## Detail of commands

**FILE < filename > -- insert lines from "filename" in current processing sequence**

**BLACKSKY -- sets entire background to black -- useful for designing new flares**

**LEAVES < N > -- sets the number of leaves on the aperture to N**

**FLARE < type radius reverse_radius [D_per_angle E_per_angle] >**

- *type* defines the flare type, 1 is simple ray, 3 gives a set of parallel rays for every leaf.
- *radius* is length of ray in pixels
- *reverse_radius* sets how long to make a the ray starting at 180 degrees from the main ray, this is a proportion which will be applied to *radius*


** (note, SCALE and ROTATE commands will overwrite previous SCALE and ROTATE values) **

**SCALE < a_h s_s s_v s_radius s_dist >**
>FOR ALL ELEMENTS (rays & ghosts) that follow in the data file:

- *a_h* is added to the hue
- element saturation is muliplied by *s_s*
- element value is multiplied by *s_v*
- element radius is multiplied by *s_radius*
- element distance is multiplied by *s_dist*


**ROTATE < D > -- rotate following elements by *D* degrees**

**IMAGE_CENTER < x y > defines the true center of the image to be x y. Units are a proportion, (default is imagewidth/2, imageheight/2)**

**SUN_CENTER < x y > defines the location of the sun (generating the rays) to be x y. Units are a proportion, (default is imagewidth/2, imageheight/2)**  The location of the sun relative to the image center changes the rendering of rays and aperture ghosts

**MOVE_CENTER < dx dy > move both IMAGE_CENTER and SUN_CENTER by dx, dy**

**GHOST_LINE < angle CA > -- creates a series of ghost-like images in a line at *angle* degrees from the center, *CA* is amount of chromatic abberation applied to the ghosts**

> additional lines that follow GHOST_LINE:

 < dist radius H S V p_poly fill >

-*dist* is the ghost distance from image center in proportion
-*radius* is the radius of the ghost in pixels
- *H,S,V*  are hue,saturation and value components of ghost color
- *p_poly* defines the shape of the ghosts, 0 to 1, 0 is circle, 1 is polygon with number of aperture leaves
- 
**GHOST_RING < n dist radius CA > -- creates a series of ghost-like images in a ring  at distance *dist* pixels  with radius *radius* pixels, *CA* is amount of chromatic abberation applied to the ghosts**

> additional lines that follow GHOST_LINE:

 < angle R G B p_poly fill >

-*angle* is the angle in degrees for the ghost from image center
- *H,S,V*  are hue,saturation and value components of ghost color
- *p_poly* defines the shape of the ghosts, 0 to 1, 0 is circle, 1 is polygon with number of aperture leaves
- *fill*  how to fill the color, 0 is a ring, 1 a disk, 4 a more gaussian like fill.  *fill* is a real number from 0 to 4

