# -s pour save la final position

from manim import *
import numpy as np
import pandas as pd
import random
from operator import itemgetter

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
        self.wait(3)


class scale_points(PCA):
    def construct(self):
        self.initial_setup()

        data = self.import_points()
        points_centered = VGroup(*[Dot([row["a_centered"], row["b_centered"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        self.add(self.plane, points_centered)
        self.play(Transform(points_centered, points_scaled, run_time = 2))
        self.wait(3)

class distance_vis(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE)

        projections = self.create_projections(points_scaled, PC1)
        
        dist_brace = Brace(Line([0, 0, 0], projections["points"][7].get_center(), color = RED), direction=DOWN)
        dist_formula = MathTex("d^2 = {:.2f}".format(self.calc_square_dist(projections["points"][7]))).next_to(dist_brace, DOWN, buff=0)


        self.initial_setup()       
        self.add(self.plane, points_scaled, PC1, projections, dist_brace, dist_formula)

class distance_vis2(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE)

        projections = self.create_projections(points_scaled, PC1)
        
        dist_brace = BraceBetweenPoints([0, 0, 0], projections["points"][7].get_center())
        dist_brace.add_updater(lambda d: d.become(BraceBetweenPoints([0, 0, 0], projections["points"][7].get_center())))
        dist_formula = MathTex("d^2 = {:.2f}".format(self.calc_square_dist(projections["points"][7]))).next_to(dist_brace, DOWN, buff=0)
        dist_formula.add_updater(lambda d: d.become(MathTex("d^2 = {:.2f}".format(self.calc_square_dist(projections["points"][7])))).next_to(dist_brace, DOWN, buff=0))


        self.initial_setup()       
        self.add(self.plane, points_scaled, PC1, projections, dist_brace, dist_formula)
        self.play(Rotate(PC1, angle = TAU/20, about_point = ORIGIN, run_time = 0.5))
        self.wait(15)


class tot_distance_calc(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])
        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE).rotate(TAU/20, about_point = ORIGIN)
        projections = self.create_projections(points_scaled, PC1)

        curr_point = projections["points"][7]

        dist_brace = BraceBetweenPoints([0, 0, 0], curr_point.get_center())
        # dist_brace.add_updater(lambda d: d.become(BraceBetweenPoints([0, 0, 0], curr_point.get_center())))
        dist_formula = MathTex("d^2_8 = ").next_to(dist_brace, DOWN, buff=0).shift(0.20 * LEFT)
        #dist_formula.add_updater(lambda d: d.become(MathTex("d^2 = ")).next_to(dist_brace, DOWN if dist_brace.get_center()[0] > 0 else UP, buff=0))

        
        dist_value = MathTex("{:.2f}".format(self.calc_square_dist(curr_point))).next_to(dist_formula, RIGHT, buff=0).shift(0.10 * RIGHT + 0.05 * DOWN)
        #dist_value.add_updater(lambda d: d.become(MathTex("{:.2f}".format(self.calc_square_dist(curr_point)))).next_to(dist_formula, RIGHT, buff=0))



        x_val = [i.get_center()[0] for i in projections["points"]]
        projections["points"] = VGroup(*[projections["points"][i] for i in [x_val.index(i) for i in sorted(x_val)]])

        cum_total = 0
        cum_formula = MathTex(r"\sum d^2_i = ").to_corner(DR).shift(1.15 * LEFT)
        cum_value = MathTex("{:.2f}".format(cum_total)).to_corner(DR).shift(0.20 * UP)


        self.initial_setup()       
        self.add(self.plane, points_scaled, PC1, projections, dist_brace, dist_formula, dist_value, cum_formula, cum_value)

        anim_run_time = 0.7

        for i in range(len(points_scaled) - 1, -1, -1):

            

            new_dist_brace = BraceBetweenPoints([0, 0, 0], projections["points"][i].get_center())
            if projections["points"][i].get_center()[0] > 0:
                new_dist_formula = MathTex("d^2_{} = ".format(i+1)).next_to(new_dist_brace, DOWN, buff=0).shift(0.20 * LEFT)
                new_dist_value = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"][i]))).next_to(new_dist_formula, RIGHT, buff=0).shift(0.10 * RIGHT + 0.05 * DOWN)
            else:
                new_dist_formula = MathTex("d^2_{} = ".format(i+1)).next_to(new_dist_brace, UP, buff=0).shift(0.20 * LEFT)
                new_dist_value = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"][i]))).next_to(new_dist_formula, RIGHT, buff=0).shift(0.10 * RIGHT + 0.05 * DOWN)

            new_dist_formula2 = MathTex("d^2_{} = ".format(i+1)).next_to(cum_formula, UP, buff=0).shift(i * 0.65 * UP + 0.2 * UP + 0.40 * RIGHT)
            new_dist_value2 = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"][i]))).next_to(new_dist_formula2, RIGHT, buff=0).shift(0.10 * RIGHT + 0.05 * DOWN)
            cum_total += self.calc_square_dist(projections["points"][i])
            new_cum_value = MathTex("{:.2f}".format(cum_total)).to_corner(DR).shift(0.20 * UP)

            trans_brace = Transform(dist_brace, new_dist_brace, run_time = anim_run_time)
            trans_formula = Transform(dist_formula, new_dist_formula, run_time = anim_run_time)
            trans_value = Transform(dist_value, new_dist_value, run_time = anim_run_time)
            create_formula = Create(new_dist_formula2, run_time = anim_run_time/100)
            create_value = Create(new_dist_value2, run_time = anim_run_time/100)
            trans_cum_value = Transform(cum_value, new_cum_value, run_time = anim_run_time)

            self.play(trans_brace, trans_formula, trans_value, create_formula, create_value, trans_cum_value)
        self.wait(4)

        
class pc1_find(PCA):
    def construct(self):        
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE)

        projections = self.create_projections(points_scaled, PC1)

        formula = MathTex(r"\sum d^2_i = ").to_corner(DR).shift(1.15 * LEFT)
        dist_indicator = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)
        dist_indicator.add_updater(lambda d: d.become(MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)))

        self.initial_setup()       
        self.add(self.plane, PC1, points_scaled, projections, formula, dist_indicator)
        self.play(Rotate(PC1, angle = TAU/2, about_point = ORIGIN, run_time = 10))
        self.wait(4)


class pc1_set(PCA):
    def construct(self):  
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE)

        projections = self.create_projections(points_scaled, PC1)

        formula = MathTex(r"\sum d^2_i = ").to_corner(DR).shift(1.15 * LEFT)
        dist_indicator = MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)
        dist_indicator.add_updater(lambda d: d.become(MathTex("{:.2f}".format(self.calc_square_dist(projections["points"]))).to_corner(DR).shift(0.20 * UP)))
        
        self.initial_setup()      
        self.add(self.plane, PC1, points_scaled, projections, formula, dist_indicator)
        self.play(Rotate(PC1, angle = TAU/8, about_point = ORIGIN, run_time = 2))
        self.wait(0.5)
        self.play(Transform(PC1, Line([-10, 0, 0], [10, 0, 0], color = PROMUTUEL_YELLOW).rotate(angle = TAU/8, about_point = ORIGIN), run_time = 0.1))
        self.wait(4)

        
class pc2_find(PCA):
    def construct(self):
        pass

class rotate_axes(PCA):
    def construct(self):
        data = self.import_points()
        points_scaled = VGroup(*[Dot([row["a_scaled"], row["b_scaled"], 0], color = PROMUTUEL_YELLOW) for _, row in data.iterrows()])

        PC1 = Line([-10, -0, 0], [10, 0, 0], color = BLUE).rotate(TAU/8, about_point = ORIGIN)
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
