# TODO: Évaluer l'utilisation de VGroup au lieu de lists pour hold des mobjects

from manim import *
import numpy as np
import pandas as pd
import random

PROMUTUEL_YELLOW = "#FFC300"
PROMUTUEL_GREY = "#adafb2"


class PCA(Scene):
    def initial_setup(self):
        self.camera.background_color = PROMUTUEL_GREY

        self.plane = NumberPlane(
            x_range=[-10, 10], 
            y_range=[-10, 10], 
            background_line_style={
                "stroke_color": PROMUTUEL_YELLOW,
                "stroke_width": 3,
                "stroke_opacity": 0.2
            }
        )
        
    def construct(self):
        pass 

    def import_points(self):
        return pd.read_csv("./data/pca_example.csv") 

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

    def create_projections(self, points, line, point_color = RED):

        # Vraiment dégueux, mais ça ne marche pas quand je loop...
        points_projected = [Dot(self.orthogonal_project_point_on_line(dot, line), color = point_color, radius = DEFAULT_DOT_RADIUS*0.75, fill_opacity=0.75) for dot in points]

        points_projected[0].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[0], line)))
        points_projected[1].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[1], line)))
        points_projected[2].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[2], line)))
        points_projected[3].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[3], line)))
        points_projected[4].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[4], line)))
        points_projected[5].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[5], line)))
        points_projected[6].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[6], line)))
        points_projected[7].add_updater(lambda m: m.move_to(self.orthogonal_project_point_on_line(points[7], line)))

        lines_projected = [DashedLine(points[i].get_center(), points_projected[i].get_center()) for i in range(8)]

        lines_projected[0].add_updater(lambda d: d.become(DashedLine(points[0].get_center(), points_projected[0].get_center())))
        lines_projected[1].add_updater(lambda d: d.become(DashedLine(points[1].get_center(), points_projected[1].get_center())))
        lines_projected[2].add_updater(lambda d: d.become(DashedLine(points[2].get_center(), points_projected[2].get_center())))
        lines_projected[3].add_updater(lambda d: d.become(DashedLine(points[3].get_center(), points_projected[3].get_center())))
        lines_projected[4].add_updater(lambda d: d.become(DashedLine(points[4].get_center(), points_projected[4].get_center())))
        lines_projected[5].add_updater(lambda d: d.become(DashedLine(points[5].get_center(), points_projected[5].get_center())))
        lines_projected[6].add_updater(lambda d: d.become(DashedLine(points[6].get_center(), points_projected[6].get_center())))
        lines_projected[7].add_updater(lambda d: d.become(DashedLine(points[7].get_center(), points_projected[7].get_center())))

        return VDict([("points", VGroup(*points_projected)), ("lines", VGroup(*lines_projected))])
class center_points(PCA):
    def construct(self):
        self.initial_setup()

        data = self.import_points()
        points_init = VGroup(*[Dot([row["a"], row["b"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        points_centered = VGroup(*[Dot([row["a_centered"], row["b_centered"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        self.add(self.plane, points_init)
        self.play(Transform(points_init, points_centered, run_time = 2))
        self.wait(15)


class scale_points(PCA):
    def construct(self):
        self.initial_setup()

        data = self.import_points()
        points_centered = VGroup(*[Dot([row["a_centered"], row["b_centered"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        self.add(self.plane, points_centered)
        self.play(Transform(points_centered, points_scaled, run_time = 2))
        self.wait(15)


class distance_vis(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        PC1 = Line([-10, -0, 0], [10, 0, 0])

        projections = self.create_projections(points_scaled, PC1)
        
        dist_brace = Brace(Line([0, 0, 0], projections["points"][2].get_center(), color = RED), direction=DOWN)
        dist_formula = MathTex("d = {:.2f}".format(self.calc_square_dist(projections["points"][2])**0.5)).next_to(dist_brace, DOWN, buff=0)


        self.initial_setup()       
        self.add(self.plane, points_scaled, PC1, projections, dist_brace, dist_formula)

class distance_vis2(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        PC1 = Line([-10, -0, 0], [10, 0, 0])

        projections = self.create_projections(points_scaled, PC1)
        
        dist_brace = BraceBetweenPoints([0, 0, 0], projections["points"][2].get_center())
        dist_brace.add_updater(lambda d: d.become(BraceBetweenPoints([0, 0, 0], projections["points"][2].get_center())))
        dist_formula = MathTex("d = {:.2f}".format(self.calc_square_dist(projections["points"][2])**(0.5))).next_to(dist_brace, DOWN, buff=0)
        dist_formula.add_updater(lambda d: d.become(MathTex("d = {:.2f}".format(self.calc_square_dist(projections["points"][2])**(0.5)))).next_to(dist_brace, DOWN, buff=0))


        self.initial_setup()       
        self.add(self.plane, points_scaled, PC1, projections, dist_brace, dist_formula)
        self.play(Rotate(PC1, angle = TAU/20, about_point = ORIGIN, run_time = 0.5))
        self.wait(60)


class pc1_find(PCA):
    def construct(self):        
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0])

        projections = self.create_projections(points_scaled, PC1)

        formula = MathTex(r"\sum d^2_i = ").to_corner(DR).shift(1.15 * LEFT)
        dist_indicator = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)
        dist_indicator.add_updater(lambda d: d.become(MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)))

        self.initial_setup()       
        self.add(self.plane, PC1, points_scaled, projections, formula, dist_indicator)
        self.wait(0.5)
        self.play(Rotate(PC1, angle = TAU/2, about_point = ORIGIN, run_time = 10))
        self.wait(0.5)
        self.play(Rotate(PC1, angle = TAU/8, about_point = ORIGIN, run_time = 2))
        self.wait(1.5)

        
class pc2_find(PCA):
    def construct(self):
        pass

class rotate_axes(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0]).rotate(TAU/8, about_point = ORIGIN)
        PC2 = Line([-10, -0, 0], [10, 0, 0]).rotate(-TAU/8, about_point = ORIGIN)

        self.initial_setup()
        
        self.add(self.plane, points_scaled, PC1, PC2)
        self.play(Uncreate(self.plane))
        self.play(Create(self.plane))

        self.play(Rotate(VGroup(*[PC1, PC2, points_scaled]), angle = -TAU/8, about_point = ORIGIN, run_time = 2))

        self.initial_setup()

        self.wait(1)
        self.play(Create(self.plane))
        self.wait(1.5)





class test(Scene):
    def construct(self):
        line = Line([-10, -0, 0], [10, 0, 0])

        line.rotate(TAU/8, about_point = ORIGIN)
        self.play(Create(line))
        self.play(Uncreate(line, remover=False))
        # line = Line([-10, -0, 0], [10, 0, 0])
        self.play(Create(line))
