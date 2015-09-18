#!/usr/bin/python

#beam, energy, run, V_Si1,V_Si2, diode type, Pb X0
CONDITIONS=[
    (11, 50,   3320, None,  -700,  "FZ 120 N",  3),
    (11, 50,   3317, None,  -600,  "FZ 120 N",  3),
    (11, 50,   3318, None,  -400,  "FZ 120 N",  3),
    (11, 50,   3319, None,  -400,  "FZ 120 N",  3),
    (11, 50,   3316, None,  -700,  "FZ 320 N",  3),
    (11, 50,   3312, None,  -600,  "FZ 320 N",  3),
    (11, 50,   3313, None,  -600,  "FZ 320 N",  3),  
    (11, 50,   3314, None,  -400,  "FZ 320 N",  3),  
    (11, 50,   3315, None,  -400,  "FZ 320 N",  3),
    (11, 50,   3363, 600,    600,  "FZ 120 P",  0),  
    (11, 50,   3323, 400,    400,  "FZ 120 P",  3),  
    (11, 50,   3321, 600,    600,  "FZ 120 P",  3),  
    (11, 50,   3322, 600,    600,  "FZ 120 P",  3),
    (11, 50,   3324, 700,    700,  "FZ 120 P",  3),  
    (11, 50,   3325, 700,    700,  "FZ 120 P",  3),
    (11, 50,   3326, 700,    700,  "FZ 120 P",  3),
    (11, 50,   3329, 700,    700,  "FZ 120 P",  3),
    (11, 50,   3358, 600,    600,  "FZ 120 P",  4),  
    (11, 50,   3362, 600,    600,  "FZ 120 P",  4),  
    (11, 50,   3369, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3370, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3371, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3372, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3373, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3374, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3375, 600,    600,  "FZ 200 P",  0),  
    (13, 150,  3376, 600,    600,  "FZ 200 P",  0),  
    (11, 50,   3296, 400,    400,  "FZ 200 P",  3),  
    (11, 50,   3295, 600,    600,  "FZ 200 P",  3),  
    (11, 50,   3365, 600,    600,  "FZ 200 P",  4),  
    (11, 50,   3367, 600,    600,  "FZ 200 P",  4),  
    (11, 50,   3346, 600,    600,  "FZ 320 P",  0),
    (11, 50,   3347, 600,    600,  "FZ 320 P",  1),  
    (11, 50,   3334, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3335, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3336, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3337, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3338, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3339, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3340, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3341, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3342, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3348, 600,    600,  "FZ 320 P",  2),  
    (11, 50,   3332, 400,    400,  "FZ 320 P",  3),  
    (11, 50,   3330, 600,    600,  "FZ 320 P",  3),  
    (11, 50,   3333, 600,    600,  "FZ 320 P",  3),  
    (11, 50,   3349, 600,    600,  "FZ 320 P",  3),  
    (11, 50,   3350, 600,    600,  "FZ 320 P",  3),
    (11, 50,   3352, 600,    600,  "FZ 320 P",  3),  
    (11, 50,   3331, 700,    700,  "FZ 320 P",  3),  
    (11, 50,   3353, 600,    600,  "FZ 320 P",  4),  
    (11, 50,   3355, 600,    600,  "FZ 320 P",  4),
    (11, 50,   3356, 600,    600,  "FZ 320 P",  4) 
]


#diode : (voltage full depletion (V_fd), thickness at V_fd, thickness at 400 V, thickness at 600V) 
DIODES={
    "FZ 320 N":(210,288,291, 292),
    "FZ 320 P":(229,282,284,285),
    "FZ 200 N":(110,208,218,220),
    "FZ 200 P":(85,200,210,211.5),
    "FZ 120 N":(100,133,145,148),
    "FZ 120 P":(71,120,131,133)
    }
