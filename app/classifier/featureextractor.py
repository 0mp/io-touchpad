from math import sqrt, atan, pi, fabs, hypot


class Point:
    """"Class represtenting a point, made to make notation more intuitive."""

    def __init__(self, x, y):
        """Create a point with two coordinates."""
        self.x_cord = x
        self.y_cord = y

    def flip_vertically(self):
        """Flip point symmetrically to OX axis."""
        self.y_cord = -self.y_cord

    def flip_horizontally(self):
        """Flip point symmetrically to OY axis."""
        self.x_cord = -self.x_cord


class Line:
    """Line is represented by two points that it connects."""

    def __init__(self, point1, point2):
        """Initialize line as connection of two points."""
        self.point1 = point1
        self.point2 = point2

    def center_point(self):
        """Calculate center of mass of a line."""
        new_point_x = float(self.point1.x_cord + self.point2.x_cord) / 2
        new_point_y = float(self.point1.y_cord + self.point2.y_cord) / 2
        return Point(new_point_x, new_point_y)

    def length(self):
        """Calculate length of the line."""
        return hypot(self.point1.x_cord - self.point2.x_cord,
                     self.point1.y_cord - self.point2.y_cord)

    def ratio_point(self, ratio):
        """Given ratio point calculates combination of points."""
        return Point(self.point1.x_cord * (1 - ratio) +
                     self.point2.x_cord * ratio,
                     self.point1.y_cord * (1 - ratio) +
                     self.point2.y_cord * ratio)


class Scaler:
    """Class used to scale points and move it to mass center."""

    def __init__(self, min_point, max_point, origin):
        """Class initiates with center of mass and scale."""
        self.min_point = Point(min_point.x_cord, min_point.y_cord)
        self.max_point = Point(max_point.x_cord, max_point.y_cord)
        self.origin = Point(origin.x_cord, origin.y_cord)

    def move_point(self, point):
        """Move point so that origin is in the center."""
        moved_point = Point(point.x_cord - self.origin.x_cord,
                            point.y_cord - self.origin.y_cord)
        return moved_point

    def scale_point(self, point):
        """Move to point to the destinated scaled place.

        Move the point to the place where it would be
        on properly scaled and moved plane
        (scale only squarely, not rectangularly).
        """
        scale = 1000
        moved_point = self.move_point(point)

        diff_x = self.max_point.x_cord - self.min_point.x_cord
        diff_y = self.max_point.y_cord - self.min_point.y_cord

        if diff_x > diff_y:
            drawn_scale = diff_x
            if diff_x != 0:
                return Point(moved_point.x_cord / drawn_scale * scale,
                             moved_point.y_cord / drawn_scale * scale)
        else:
            drawn_scale = diff_y
            if diff_y != 0:
                return Point(moved_point.x_cord / drawn_scale * scale,
                             moved_point.y_cord / drawn_scale * scale)
        return Point(0, 0)


class Curve:
    """Represent only curves made of linear segments."""

    def __init__(self, starting_point):
        """Initiate curve with beginning point."""
        self.length = 0
        self.list_of_points = []
        self.colors = []

        self.list_of_points.append(starting_point)
        self.number_of_points = 1  # necessary to calculate new center of mass
        self.center_of_mass = \
            Point(starting_point.x_cord, starting_point.y_cord)

    def actualise_center_of_mass(self, point, line_length):
        """Actualize center of mass acording to point and length."""
        if self.length + line_length > 0:

            new_line_mass_x = point.x_cord * line_length
            whole_mass_x = self.center_of_mass.x_cord * self.length
            self.center_of_mass.x_cord = \
                (whole_mass_x + new_line_mass_x) / \
                (self.length + line_length)

            new_line_mass_y = point.y_cord * line_length
            whole_mass_y = self.center_of_mass.y_cord * self.length
            self.center_of_mass.y_cord = \
                (whole_mass_y + new_line_mass_y) / \
                (self.length + line_length)

    def add_point(self, point):
        """Add point to curve taking care of actualising parameters."""
        last_point = self.list_of_points[len(self.list_of_points) - 1]

        self.list_of_points.append(point)
        added_line = Line(last_point, point)
        line_length = added_line.length()
        new_line_center_point = added_line.center_point()

        self.actualise_center_of_mass(new_line_center_point, line_length)
        self.length += line_length

    def hard_add_point(self, point):
        """Add point without actualising parameters."""
        self.list_of_points.append(point)

    def attach_colors(self, colors):
        """Attache list of colors to the curve."""
        self.colors = colors

    def add_color(self, color):
        """Add color to the list of point of the curve."""
        self.colors.append(color)


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
            min_point.x_cord = point.x_cord
            min_point.y_cord = point.y_cord
            max_point.x_cord = point.x_cord
            max_point.y_cord = point.y_cord
            initiated_starting_values = True

        if point.x_cord < min_point.x_cord:
            min_point.x_cord = point.x_cord
        if point.x_cord > max_point.x_cord:
            max_point.x_cord = point.x_cord
        if point.y_cord < min_point.y_cord:
            min_point.y_cord = point.y_cord
        if point.y_cord > max_point.y_cord:
            max_point.y_cord = point.y_cord

    return min_point, max_point


def create_curve(signal_list):
    if not signal_list:
        return
    curve = Curve(Point(signal_list[0].get_x(), signal_list[0].get_y()))

    for i in range(1, len(signal_list)):
        point = Point(signal_list[i].get_x(), signal_list[i].get_y())
        curve.add_point(point)

    return curve


def create_normalized_curve(curve, min_point, max_point, colors):
    # creates list of equdistant NUMBER_OF_POINTS points
    # that represents the same shape as signal_list list

    length_of_one_line = curve.length / (NUMBER_OF_POINTS - 1)
    # there is one more point then the number of lines
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
        while travelled_distance > length_of_one_line:
            # there should be a new points between these two
            travelled_distance = (travelled_distance - length_of_one_line)
            overdue = section.length() - travelled_distance
            # section_length = overdue + x
            # overdue is included in previous line,
            # so x must be added to the new line distance

            # determining new point coordinates
            ratio = overdue / section.length()
            point = section.ratio_point(ratio)
            scaled_point = scaler.scale_point(point)
            section = Line(point, next_point)
            # in case there should be more points added then one here

            normalized_curve.hard_add_point(scaled_point)

            normalized_curve.add_color(colors[i])

    while len(normalized_curve.list_of_points) < NUMBER_OF_POINTS:
        normalized_curve.hard_add_point(
            curve.list_of_points[len(curve.list_of_points) - 1])
        normalized_curve.add_color(colors[len(curve.list_of_points) - 1])

    normalized_curve.center_of_mass = Point(0, 0)
    normalized_curve.length = curve.length
    return normalized_curve


def draw_new_points(list_of_points):
    # testing function, to use with matrixanalyser
    list_of_signal_points, colors = filter_points_from_signals(list_of_points)
    new_curve = normalize_points(list_of_signal_points, colors)
    new_points = new_curve.list_of_points
    with open('tools/matrixanalyser/data/coordinates2.data', 'w') \
            as drawing_file:
        for point in new_points:
            drawing_file.write("%d %d\n" % (point[0], point[1]))


def angle_between_line_and_xaxis(point1, point2):
    """Calculate angle beetwen xaxis and line.

    Xaxis is joint to point2, angle on the left side.
    """
    if point2.x_cord != point1.x_cord:
        return atan((point2.y_cord - point1.y_cord) /
                    (point2.x_cord - point1.x_cord))
    if point2.y != point1.y:
        return (pi / 2) * (point2.y_cord - point1.y_cord) / \
            abs(point2.y_cord - point1.y_cord)
    return 0


def dot_product(vector1, vector2):
    return vector1[0] * vector2[0] + vector1[1] + vector2[1]


def get_angle_list(curve):
    # gets list of angles between lines and x_axis
    feature_list = []
    list_of_points = curve.list_of_points

    for i in range(len(list_of_points) - 1):
        point = list_of_points[i]
        next_point = list_of_points[i + 1]
        angle = angle_between_line_and_xaxis(point, next_point)

        # scaling angle
        angle = 2 * angle / pi * (SCALE / ANGLE_DOWNSCALE)

        # appends absolute value of the angle
        feature_list.append(fabs(angle))
    return feature_list


def join_features(list_of_points, list_of_feature1, colors):
    """Function used to join features into single list.

    Assumes length of points is the biggest here
    # join lists of coordinates with features (one feature for now),
    # in order x,y,feature1,x,y,feature1,....
    """
    feature_list = []
    feature1_length = len(list_of_feature1)
    for i in range(len(list_of_points)):
        point = list_of_points[i]
        feature_list.append(point.x_cord)
        feature_list.append(point.y_cord)
        feature_list.append(colors[i])
        if i < feature1_length:
            feature_list.append(list_of_feature1[i])
    return feature_list


def normalize_points(list_of_signal_points, colors):
    # returns list of points that represents the same
    # figure but is in our preferred standard
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
            while j < length and not \
                    list_of_signals[j].is_raising_finger_signal()\
                    and not list_of_signals[j].is_proper_signal_of_point():
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
