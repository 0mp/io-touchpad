from math import pow, sqrt, acos, atan, pi, fabs


def center_of_line(point1, point2):
    return (point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2


def length_of_line(point1, point2):
    return sqrt(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2))


SCALE = 1000
NUMBER_OF_POINTS = 40
ANGLE_DOWNSCALE = 28
COLOR_DOWNSCALE = 2


def calculate_border_points(signal_list):
    initiated_starting_values = False;
    for signal in signal_list:
        point = signal.get_x(), signal.get_y()
        if not initiated_starting_values:
            min_x = point[0]
            max_x = point[0]
            min_y = point[1]
            max_y = point[1]
            initiated_starting_values = True

        if point[0] < min_x:
            min_x = point[0]
        if point[0] > max_x:
            max_x = point[0]
        if point[1] < min_y:
            min_y = point[1]
        if point[1] > max_y:
            max_y = point[1]
    return min_x, min_y, max_x, max_y


def calculate_center_of_mass_and_length(signal_list):
    curve_length = 0
    whole_mass_x = 0
    whole_mass_y = 0
    for i in range(len(signal_list) - 1):
        point = signal_list[i].get_x(), signal_list[i].get_y()
        next_point = signal_list[i + 1].get_x(), signal_list[i + 1].get_y()

        # curve length
        length = length_of_line(point, next_point)
        curve_length += length

        # adding to center of mass parameters
        center_x, center_y = center_of_line(point, next_point)
        whole_mass_x += center_x * length
        whole_mass_y += center_y * length

    if curve_length == 0: return (signal_list[0].get_x(), signal_list[0].get_y()), 0

    center_of_mass = whole_mass_x / curve_length, whole_mass_y / curve_length
    return center_of_mass, curve_length


def ratio_point_of_line(point1, point2, ratio):
    return point1[0] * (1 - ratio) + point2[0] * ratio, point1[1] * (1 - ratio) + point2[1] * ratio


def scale_point(point, min_x, min_y, max_x, max_y, origin):
    # moves the point to the place where it would be on properly scaled and moved plane
    # (scales only squarely, not rectangularly)
    moved_point = point[0] - origin[0], point[1] - origin[1]
    diff_x = max_x - min_x
    diff_y = max_y - min_y
    if diff_x > diff_y:
        drawn_scale = diff_x
        if diff_x != 0:
            return moved_point[0] / drawn_scale * SCALE, moved_point[1] / drawn_scale * SCALE
    else:
        drawn_scale = diff_y
        if diff_y != 0:
            return moved_point[0] / drawn_scale * SCALE, moved_point[1] / drawn_scale * SCALE
    return 0, 0


def create_normalized_list_of_points(center_of_mass, min_x, min_y, max_x, max_y, curve_length, signal_list, colors):
    # creates list of equdistant NUMBER_OF_POINTS points that represents the same shape as signal_list list

    length_of_one_line = (curve_length) / NUMBER_OF_POINTS
    # curve_length-1 to be sure there are NUMBER_OF_POINTS points

    new_points = []
    new_colors = []

    travelled_distance = 0
    for i in range(len(signal_list) - 1):
        point = signal_list[i].get_x(), signal_list[i].get_y()
        next_point = signal_list[i + 1].get_x(), signal_list[i + 1].get_y()
        section_length = length_of_line(point, next_point)
        travelled_distance += section_length
        while travelled_distance > length_of_one_line:  # there should be a new points between these two
            travelled_distance = (travelled_distance - length_of_one_line)
            overdue = section_length - travelled_distance
            # section_length = overdue + x
            # overdue is included in previous line, so x must be added to the new line distance

            # determining new point coordinates
            ratio = overdue / section_length
            point = ratio_point_of_line(point, next_point,
                                        ratio)  # in case there should be more points added then one here
            section_length = length_of_line(point, next_point)
            new_points.append(scale_point(point, min_x, min_y, max_x, max_y, center_of_mass))
            new_colors.append(colors[i])
            # new_points.append(point)

    while len(new_points) != NUMBER_OF_POINTS:
        point = signal_list[len(signal_list) - 1].get_x(), signal_list[len(signal_list) - 1].get_y()
        new_points.append(scale_point(point, min_x, min_y, max_x, max_y, center_of_mass))
        new_colors.append(colors[len(signal_list) - 1])
    return new_points, new_colors


def draw_new_points(list_of_points):
    #testing function, to use with matrixanalyser
    center_of_mass,minX,minY,maxX,maxY,curve_length = calculate_normalization_values(list_of_points)
    new_points, _ = create_normalized_list_of_points(center_of_mass,minX,minY,maxX,maxY,curve_length,list_of_points)

    get_features(list_of_points)

    with open('tools/matrixanalyser/data/coordinates2.data', 'w') as drawing_file:
        for point in new_points:
            drawing_file.write("%d %d\n" % (point[0], point[1]))


def get_angle_between_line_and_xaxis(point1, point2):  # xaxis joint to point2, angle on the left side
    if point2[0] != point1[0]:
        return atan((point2[1] - point1[1]) / (point2[0] - point1[0]))
    if point2[1] != point1[1]:
        return (pi / 2) * (point2[1] - point1[1]) / abs(point2[1] - point1[1])
    return 0


def dot_product(vector1, vector2):
    return vector1[0] * vector2[0] + vector1[1] + vector2[1]


def angle_between_lines(point1, point2, point3):  # point2 is middle point
    vector1 = point1[0] - point2[0], point1[1] - point2[1]
    vector2 = point3[0] - point2[0], point3[1] - point2[1]
    # todo: wykrywanie czy nie trzeba zamienić wyniku na 2pi - tenCoTerazWychodzi
    return acos(dot_product(vector1, vector2) / (length_of_line(point1, point2) * length_of_line(point2, point3)))


def get_angle_list(list_of_points):
    # gets list of angles between lines and x_axis
    feature_list = []

    for i in range(len(list_of_points) - 1):
        point = list_of_points[i]
        next_point = list_of_points[i + 1]
        angle = get_angle_between_line_and_xaxis(point, next_point)

        # scaling angle
        angle = 2 * angle / pi * (SCALE / ANGLE_DOWNSCALE)

        # appends absolute value of the angle
        feature_list.append(fabs(angle))
    return feature_list


def join_features(list_of_points, list_of_feature1, colors):  # assumes length of points is the biggest here
    # join lists of coordinates with features (one feature for now), in order x,y,feature1,x,y,feature1,....
    feature_list = []
    feature1_length = len(list_of_feature1)
    for i in range(len(list_of_points)):
        point = list_of_points[i]
        feature_list.append(point[0])
        feature_list.append(point[1])
        feature_list.append(colors[i])
        if i < feature1_length:
            feature_list.append(list_of_feature1[i])
    return feature_list


def normalize_points(list_of_points, colors):
    # returns list of points that represents the same figure but is in our preferred standard
    min_x, min_y, max_x, max_y = calculate_border_points(list_of_points)
    center_of_mass, curve_length = calculate_center_of_mass_and_length(list_of_points)

    new_points, new_colors = create_normalized_list_of_points(center_of_mass,minX,minY,maxX,maxY,curve_length,list_of_points,colors)
    return new_points, new_colors


def filter_points_from_signals(list_of_signals):
    points = []
    colors = []
    length = len(list_of_signals)
    color_scaled = SCALE / COLOR_DOWNSCALE
    for i in range(length):
        touchpad_signal = list_of_signals[i]
        if touchpad_signal.is_proper_signal_of_point():
            points.append(touchpad_signal)
            j = i + 1
            while j < length and not list_of_signals[j].is_raising_finger_signal() and \
                    not list_of_signals[j].is_proper_signal_of_point():
                j += 1
            if j == length or list_of_signals[j].is_proper_signal_of_point():
                colors.append(0)
            else:
                colors.append(color_scaled)
    return points, colors


def get_features(list_of_signals):
    # returns list of features for list of points taken from evtest
    # changing format of points to preferred
    list_of_points, colors = filter_points_from_signals(list_of_signals)
    new_points, new_colors = normalize_points(list_of_points, colors)
    # getting features different then coordinates
    angles = get_angle_list(new_points)

    # joining all the features together
    feature_list = join_features(new_points, angles, new_colors)
    return feature_list
