from math import pow,sqrt,acos,atan,pi,fabs


def center_of_line(point1,point2):
    return (point1[0]+point2[0])/2,(point1[1]+point2[1])/2


def length_of_line(point1,point2):
    return sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2))


SCALE = 1000
NUMBER_OF_POINTS = 40
ANGLE_DOWNSCALE = 20


def calculate_normalization_values(signal_list):
   #calculating values necessary for normalization (center of mass, boundary and curve length for now)

   curve_length = 0
   whole_mass_x = 0
   whole_mass_y = 0
   for i in range(len(signal_list)-1):
        point = signal_list[i].get_x(),signal_list[i].get_y()
        next_point = signal_list[i+1].get_x(),signal_list[i+1].get_y()
        if (i==0):
            minX=point[0]
            maxX=point[0]
            minY=point[1]
            maxY=point[1]
        #curve length
        length=length_of_line(point,next_point)
        curve_length += length

        #max and min parameters
        if point[0]<minX:
            minX = point[0]
        if point[0]>maxX:
            maxX = point[0]
        if point[1]<minY:
            minY = point[1]
        if point[1]>maxY:
            maxY = point[1]

        #adding to center of mass parameters
        center_x,center_y=center_of_line(point,next_point)
        whole_mass_x += center_x*length
        whole_mass_y += center_y*length


   #last point was not included yet
   point= signal_list[len(signal_list)-1].get_x(),signal_list[len(signal_list)-1].get_y() #last point can change min and max values
   if(len(signal_list)==1):
       minX = point[0]
       maxX = point[0]
       minY = point[1]
       maxY = point[1]
   if point[0]<minX:
       minX = point[0]
   if point[0]>maxX:
       maxX = point[0]
   if point[1]<minY:
       minY = point[1]
   if point[1]>maxY:
       maxY = point[1]

   if (len(signal_list)==1): return (minX,minY),minX,minY,maxX,maxY,0
   center_of_mass = whole_mass_x/curve_length,whole_mass_y/curve_length
   return center_of_mass,minX,minY,maxX,maxY,curve_length


def ratio_point_of_line(point1,point2,ratio):
    return point1[0]*(1-ratio) + point2[0]*ratio,point1[1]*(1-ratio)+ point2[1]*ratio


def scale_point(point,minX,minY,maxX,maxY,origin):
    #moves the point to the place where it would be on properly scaled and moved plane (scales only squarely, not rectangularly)
    moved_point = point[0] - origin[0],point[1] - origin[1]
    diffX = maxX - minX
    diffY = maxY - minY
    if(diffX>diffY):
        drawn_scale = diffX
        if(diffX!=0):
            return moved_point[0]/drawn_scale*SCALE,moved_point[1]/drawn_scale*SCALE
    else:
        drawn_scale = diffY
        if(diffY!=0):
            return moved_point[0]/drawn_scale*SCALE,moved_point[1]/drawn_scale*SCALE
    return 0,0

def create_normalized_list_of_points(center_of_mass, minX, minY, maxX, maxY, curve_length, signal_list): #list of signals from evtest
    #creates list of equdistant NUMBER_OF_POINTS points that represents the same shape as signal_list list

    length_of_one_line = (curve_length)/NUMBER_OF_POINTS #curve_length-1 to be sure there are NUMBER_OF_POINTS points
    new_points = []

    travelled_distance = 0
    for i in range(len(signal_list)-1):
        point = signal_list[i].get_x(),signal_list[i].get_y()
        next_point = signal_list[i+1].get_x(),signal_list[i+1].get_y()
        section_length = length_of_line(point,next_point)
        travelled_distance += section_length
        while travelled_distance > length_of_one_line: #there should be a new points between these two
            travelled_distance=( travelled_distance -length_of_one_line)
            overdue = section_length - travelled_distance
            # section_length = overdue + x
            # overdue is included in previous line, so x must be added to the new line distance

            #determining new point coordinates
            ratio = overdue/section_length
            point = ratio_point_of_line(point,next_point,ratio) #in case there should be more points added then one here
            section_length = length_of_line(point,next_point)
            new_points.append(scale_point(point,minX,minY,maxX,maxY,center_of_mass))
            #new_points.append(point)

    while(len(new_points)!=NUMBER_OF_POINTS):
        point = signal_list[len(signal_list)-1].get_x(),signal_list[len(signal_list)-1].get_y()
        new_points.append(scale_point(point,minX,minY,maxX,maxY,center_of_mass))
    return new_points

def draw_new_points(list_of_points):
    #testing function, to use with matrixanalyser
    center_of_mass,minX,minY,maxX,maxY,curve_length = calculate_normalization_values(list_of_points)
    new_points = create_normalized_list_of_points(center_of_mass,minX,minY,maxX,maxY,curve_length,list_of_points)

    get_features(list_of_points)

    with open('tools/matrixanalyser/data/coordinates2.data','w') as file:
        for point in new_points:
            file.write("%d %d\n" % (point[0], point[1]))


def get_angle_between_line_and_xaxis(point1,point2): #xaxis joint to point2, angle on the left side
    if((point2[0]!=point1[0])):
        return atan((point2[1]-point1[1])/(point2[0]-point1[0]))
    if((point2[1]!=point1[1])):
        return (pi/2)*(point2[1]-point1[1])/abs(point2[1]-point1[1])
    return 0

def dotProduct(vector1,vector2):
    return vector1[0]*vector2[0]+vector1[1]+vector2[1]


def angle_between_lines(point1,point2,point3): #point2 is middle point
    vector1 = point1[0]-point2[0],point1[1]-point2[1]
    vector2 = point3[0]-point2[0],point3[1]-point2[1]
    #todo: wykrywanie czy nie trzeba zamienić wyniku na 2pi - tenCoTerazWychodzi
    return acos(dotProduct(vector1,vector2)/(length_of_line(point1,point2)*length_of_line(point2,point3)))


def get_angle_list(list_of_points):
    #gets list of angles between lines and x_axis
    feature_list = []

    for i in range(len(list_of_points)-1):
         point= list_of_points[i]
         next_point = list_of_points[i+1]
         angle = get_angle_between_line_and_xaxis(point,next_point)

         #scaling angle
         angle = angle/pi*(SCALE/ANGLE_DOWNSCALE)

         #appends absolute value of the angle
         feature_list.append(fabs(angle))
    return feature_list


def join_features(list_of_points,list_of_feature1): #assumes length of points is the biggest here
    #join lists of coordinates with features (one feature for now), in order x,y,feature1,x,y,feature1,....
    feature_list = []
    feature1_length = len(list_of_feature1)
    for i in range(len(list_of_points)):
        point = list_of_points[i]
        feature_list.append(point[0])
        feature_list.append(point[1])
        if (i <feature1_length):
            feature_list.append(list_of_feature1[i])
    return feature_list


def get_features(list_of_points):
    #returns list of features for list of points taken from evtest
    center_of_mass,minX,minY,maxX,maxY,curve_length = calculate_normalization_values(list_of_points)
    new_points = create_normalized_list_of_points(center_of_mass,minX,minY,maxX,maxY,curve_length,list_of_points)
    angles = get_angle_list(new_points)
    feature_list = join_features(new_points,angles)
    print(len(feature_list))

    return feature_list
