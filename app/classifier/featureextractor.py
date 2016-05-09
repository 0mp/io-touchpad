from math import pow, sqrt, atan, pi, fabs, hypot


class Point:
    """" class represtenting a point, made to make notation more intuitive """

    def __init__(self, x, y):
        """ creates a point with two coordinates """
        self.x = x
        self.y = y

    def flip_vertically(self):
        """ flips point symmetrically to OX axis """
        self.y = -self.y

    def flip_horizontally(self):
        """ flips point symmetrically to OY axis """
        self.x = -self.x


class Line:
    # line is represented by two points that it connects

    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2

    def center_point(self):
        new_point_x = float(self.point1.x + self.point2.x) / 2
        new_point_y = float(self.point1.y + self.point2.y) / 2
        return Point(new_point_x, new_point_y)

    def length(self):
        return hypot(self.point1.x - self.point2.x, self.point1.y - self.point2.y)

    def ratio_point(self, ratio):
        return Point(self.point1.x * (1 - ratio) + self.point2.x * ratio,
                     self.point1.y * (1 - ratio) + self.point2.y * ratio)


class Scaler:
    def __init__(self, min_point, max_point, origin):
        self.min_point = Point(min_point.x, min_point.y)
        self.max_point = Point(max_point.x, max_point.y)
        self.origin = Point(origin.x, origin.y)

    def scale_point(self, point):
        # moves the point to the place where it would be on properly scaled and moved plane
        # (scales only squarely, not rectangularly)
        scale = 1000
        moved_point = Point(point.x - self.origin.x, point.y - self.origin.y)
        diff_x = self.max_point.x - self.min_point.x
        diff_y = self.max_point.y - self.min_point.y

        if diff_x > diff_y:
            drawn_scale = diff_x
            if diff_x != 0:
                return Point(moved_point.x / drawn_scale * scale, moved_point.y / drawn_scale * scale)
        else:
            drawn_scale = diff_y
            if diff_y != 0:
                return Point(moved_point.x / drawn_scale * scale, moved_point.y / drawn_scale * scale)
        return Point(0, 0)


class Curve:
    # represents only curves made of linear segments

    def __init__(self, starting_point):
        self.length = 0
        self.list_of_points = []
        self.colors = []

        self.list_of_points.append(starting_point)
        self.number_of_points = 1  # necessary to calculate new center of mass
        self.center_of_mass = Point(starting_point.x, starting_point.y)

    def actualise_center_of_mass(self, point, line_length):
        if self.length + line_length > 0:
            self.center_of_mass.x = (self.center_of_mass.x * self.length + point.x * line_length) / \
                                    (self.length + line_length)
            self.center_of_mass.y = (self.center_of_mass.y * self.length + point.y * line_length) / \
                                    (self.length + line_length)

    def add_point(self, point):
        last_point = self.list_of_points[len(self.list_of_points) - 1]

        self.list_of_points.append(point)
        added_line = Line(last_point, point)
        line_length = added_line.length()
        new_line_center_point = added_line.center_point()

        self.actualise_center_of_mass(new_line_center_point, line_length)
        self.length += line_length

    def hard_add_point(self, point):
        self.list_of_points.append(point)

    def attach_colors(self, colors):
        self.colors = colors

    def add_color(self, color):
        self.colors.append(color)


class SymbolFeatures:
    # collects normalized info about symbol, ready to compare with others

    def __init__(self, curve):
        self.points = curve


def center_of_line(point1, point2):
    return (point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2


def length_of_line(point1, point2):
    return sqrt(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2))


SCALE = 1000
NUMBER_OF_POINTS = 40
ANGLE_DOWNSCALE = 28
COLOR_DOWNSCALE = 2


def calculate_border_points(signal_list):
    initiated_starting_values = False
    max_point = Point(0, 0)
    min_point = Point(0, 0)
    for signal in signal_list:
        point = Point(signal.get_x(), signal.get_y())

        if not initiated_starting_values:
            min_point.x = point.x
            min_point.y = point.y
            max_point.x = point.x
            max_point.y = point.y
            initiated_starting_values = True

        if point.x < min_point.x:
            min_point.x = point.x
        if point.x > max_point.x:
            max_point.x = point.x
        if point.y < min_point.y:
            min_point.y = point.y
        if point.y > max_point.y:
            max_point.y = point.y

    return min_point, max_point


def create_curve(signal_list):
    if not signal_list:
        return
    curve = Curve(Point(signal_list[0].get_x(), signal_list[0].get_y()))

    for i in range(1, len(signal_list)):
        point = Point(signal_list[i].get_x(), signal_list[i].get_y())
        curve.add_point(point)

    return curve


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

    if curve_length == 0:
        return (signal_list[0].get_x(), signal_list[0].get_y()), 0

    center_of_mass = whole_mass_x / curve_length, whole_mass_y / curve_length
    return center_of_mass, curve_length


def ratio_point_of_line(point1, point2, ratio):
    return point1[0] * (1 - ratio) + point2[0] * ratio, \
           point1[1] * (1 - ratio) + point2[1] * ratio


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


def create_normalized_curve(curve, min_point, max_point, colors):
    # creates list of equdistant NUMBER_OF_POINTS points that represents the same shape as signal_list list

    length_of_one_line = curve.length / (NUMBER_OF_POINTS - 1)  # there is one more point then the number of lines
    # curve_length-1 to be sure there are NUMBER_OF_POINTS points

    travelled_distance = 0
    scaler = Scaler(min_point, max_point, curve.center_of_mass)
    point = curve.list_of_points[0]

    normalized_curve = Curve(scaler.scale_point(point))
    normalized_curve.add_color(colors[0])

    for i in range(0, len(curve.list_of_points) - 1):
        point = curve.list_of_points[i]
        next_point = curve.list_of_points[i + 1]
        section = Line(point, next_point)

        travelled_distance += section.length()
        while travelled_distance > length_of_one_line:  # there should be a new points between these two
            travelled_distance = (travelled_distance - length_of_one_line)
            overdue = section.length() - travelled_distance
            # section_length = overdue + x
            # overdue is included in previous line, so x must be added to the new line distance

            # determining new point coordinates
            ratio = overdue / section.length()
            point = section.ratio_point(ratio)
            scaled_point = scaler.scale_point(point)
            section = Line(point, next_point)  # in case there should be more points added then one here
            normalized_curve.hard_add_point(scaled_point)

            normalized_curve.add_color(colors[i])

    while len(normalized_curve.list_of_points) < NUMBER_OF_POINTS:
        normalized_curve.hard_add_point(curve.list_of_points[len(curve.list_of_points) - 1])
        normalized_curve.add_color(colors[len(normalized_curve.list_of_points)])

    normalized_curve.center_of_mass = Point(0, 0)
    normalized_curve.length = curve.length
    return normalized_curve


def draw_new_points(list_of_points):
    # testing function, to use with matrixanalyser
    list_of_signal_points, colors = filter_points_from_signals(list_of_points)
    new_curve = normalize_points(list_of_signal_points, colors)
    new_points = new_curve.list_of_points
    with open('tools/matrixanalyser/data/coordinates2.data', 'w') as drawing_file:
        for point in new_points:
            drawing_file.write("%d %d\n" % (point[0], point[1]))


def get_angle_between_line_and_xaxis(point1, point2):  # xaxis joint to point2, angle on the left side
    if point2.x != point1.x:
        return atan((point2.y - point1.y) / (point2.x - point1.x))
    if point2.y != point1.y:
        return (pi / 2) * (point2.y - point1.y) / abs(point2.y - point1.y)
    return 0


def dot_product(vector1, vector2):
    return vector1[0] * vector2[0] + vector1[1] + vector2[1]


"""""
def angle_between_lines(point1, point2, point3):  # point2 is middle point
    vector1 = point1[0] - point2[0], point1[1] - point2[1]
    vector2 = point3[0] - point2[0], point3[1] - point2[1]
    return acos(dot_product(vector1, vector2) /
                (length_of_line(point1, point2) * length_of_line(point2, point3)))
"""""


def get_angle_list(curve):
    # gets list of angles between lines and x_axis
    feature_list = []
    list_of_points = curve.list_of_points

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
        feature_list.append(point.x)
        feature_list.append(point.y)
        feature_list.append(colors[i])
        if i < feature1_length:
            feature_list.append(list_of_feature1[i])
    return feature_list


def normalize_points(list_of_signal_points, colors):
    # returns list of points that represents the same figure but is in our preferred standard
    min_point, max_point = calculate_border_points(list_of_signal_points)

    curve = create_curve(list_of_signal_points)

    new_curve = create_normalized_curve(curve, min_point, max_point, colors)
    return new_curve


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


def get_new_points(list_of_signals):
    list_of_signal_points, colors = filter_points_from_signals(list_of_signals)
    new_curve = normalize_points(list_of_signal_points, colors)

    # temporary todo:zmienić
    new_points = new_curve.list_of_points
    return new_points


def get_features(list_of_signals):
    # returns list of features for list of points taken from evtest
    # changing format of points to preferred
    list_of_signal_points, colors = filter_points_from_signals(list_of_signals)

    new_curve = normalize_points(list_of_signal_points, colors)

    # temporary
    # todo: change it to objet oriented
    new_points = new_curve.list_of_points
    new_colors = new_curve.colors

    # getting features different then coordinates
    angles = get_angle_list(new_curve)

    # joining all the features together
    feature_list = join_features(new_points, angles, new_colors)
    return feature_list


def test_curve_class():
    point = Point(0, 0)
    curve = Curve(point)
    for i in range(1, 40):
        curve.add_point(Point(i, i))
    return

#
# test_curve_class();
