# TODO: Évaluer l'utilisation de VGroup au lieu de lists pour hold des mobjects

from manim import *
import numpy as np
import pandas as pd
import random

PROMUTUEL_YELLOW = "#FFC300"
PROMUTUEL_GREY = "#7F7F7F"


class PCA(Scene):
    def initial_setup(self):
        self.camera.background_color = PROMUTUEL_GREY

        plane = NumberPlane(
            x_range=[-10, 10], 
            y_range=[-10, 10], 
            background_line_style={
                "stroke_color": PROMUTUEL_YELLOW,
                "stroke_width": 3,
                "stroke_opacity": 0.2
            }
        )

        self.add(plane)
    def construct(self):
        pass 

    def import_points(self):
        return pd.read_csv("./data/pca_example.csv") 

    def generate_random_dot(self, lbound = 0, ubound = 1, generate_z = True):
        return Dot([random.uniform(lbound, ubound), random.uniform(lbound, ubound), generate_z * random.uniform(lbound, ubound)])

    def generate_random_dots(self, n, lbound = 0, ubound = 1, generate_z = True):
        return [self.generate_random_dot(lbound, ubound, generate_z) for i in range(n)]

    def create_transformation_list(self, from_list, to_list):
        return [Transform(from_list[i], to_list[i]) for i in range(len(from_list))]

    def apply_to_all(self, function, list, *args, **kwargs):
        return [function(list[i], *args, **kwargs) for i in range(len(list))]

    def orthogonal_project_point_on_line(self, dot, line):
        """Returns the point on the line segment ab that is closest to p."""
        ap = dot.get_center() - line.get_start()
        ab = line.get_end() - line.get_start()
        return line.get_start() + np.dot(ap, ab) / np.dot(ab, ab) * ab

    def calc_square_dist(self, dots):
        return sum([sum(dot.get_center()**2) for dot in dots])

class center_points(PCA):
    def construct(self):
        self.initial_setup()

        data = self.import_points()
        points_init = VGroup(*[Dot([row["a"], row["b"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        points_centered = VGroup(*[Dot([row["a_centered"], row["b_centered"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        self.add(points_init)
        self.play(Transform(points_init, points_centered, run_time = 2))

class scale_points(PCA):
    def construct(self):
        self.initial_setup()

        data = self.import_points()
        points_centered = VGroup(*[Dot([row["a_centered"], row["b_centered"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        self.add(points_centered)
        self.play(Transform(points_centered, points_scaled, run_time = 2))


class pc1_find(PCA):
    def construct(self):        
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0])

        points_projected = VGroup(*[Dot(self.orthogonal_project_point_on_line(dot, PC1), color = RED, radius = DEFAULT_DOT_RADIUS*0.75, fill_opacity=0.75) for dot in points_scaled])
        lines_projected = VGroup(*[DashedLine(point, proj_point, fill_opacity = 0.75) for point, proj_point in zip(points_scaled, points_projected)])

        # Gross, mais ça marche pas quand je loop and I don't care to fix it
        points_projected[0].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[0], PC1)))
        points_projected[1].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[1], PC1)))
        points_projected[2].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[2], PC1)))
        points_projected[3].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[3], PC1)))
        points_projected[4].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[4], PC1)))
        points_projected[5].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[5], PC1)))
        points_projected[6].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[6], PC1)))
        points_projected[7].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points_scaled[7], PC1)))

        lines_projected[0].add_updater(lambda d: d.become(DashedLine(points_scaled[0].get_center(), points_projected[0].get_center())))
        lines_projected[1].add_updater(lambda d: d.become(DashedLine(points_scaled[1].get_center(), points_projected[1].get_center())))
        lines_projected[2].add_updater(lambda d: d.become(DashedLine(points_scaled[2].get_center(), points_projected[2].get_center())))
        lines_projected[3].add_updater(lambda d: d.become(DashedLine(points_scaled[3].get_center(), points_projected[3].get_center())))
        lines_projected[4].add_updater(lambda d: d.become(DashedLine(points_scaled[4].get_center(), points_projected[4].get_center())))
        lines_projected[5].add_updater(lambda d: d.become(DashedLine(points_scaled[5].get_center(), points_projected[5].get_center())))
        lines_projected[6].add_updater(lambda d: d.become(DashedLine(points_scaled[6].get_center(), points_projected[6].get_center())))
        lines_projected[7].add_updater(lambda d: d.become(DashedLine(points_scaled[7].get_center(), points_projected[7].get_center())))

        formula = MathTex(r"\sum d^2_i = ").to_corner(DR).shift(1.15 * LEFT)
        dist_indicator = MathTex("{:.2f}".format(self.calc_square_dist(points_projected))).to_corner(DR).shift(0.20 * UP)
        #dist_indicator.add_updater(lambda d: d.set_tex_string("$\\sum d^2_i = " + str(round(self.calc_square_dist(points_projected), 2)) + "$"))
        dist_indicator.add_updater(lambda d: d.become(MathTex("{:.2f}".format(self.calc_square_dist(points_projected))).to_corner(DR).shift(0.20 * UP)))

        self.initial_setup()       
        self.add(PC1, points_scaled, points_projected, lines_projected, formula, dist_indicator)
        self.wait(0.5)
        self.play(Rotate(PC1, angle = TAU/2, about_point = ORIGIN, run_time = 10))
        

