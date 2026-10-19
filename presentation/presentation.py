from pathlib import Path

import numpy as np
from manim import *
from manim_slides import Slide, ThreeDSlide

# 4:3 frame: x in [-4, 4], y in [-3, 3]. Keep the pixel height from the
# quality flag (-ql/-qm/-qh/-qk) and derive the width from it.
config.frame_width = 8
config.frame_height = 6
config.pixel_width = config.pixel_height * 4 // 3

ASSETS = Path(__file__).parent / "assets"


def elliptic_curve(x_max=config.frame_width / 2 + 0.2, height=None, t_step=0.01, **kwargs):
    """The curve 7y^2 = (x-1)^3 + 3(x-1)^2 + 6 for x <= x_max, traced from the
    lower branch to the upper branch.

    With u = x - 1 and r the real root of u^3 + 3u^2 + 6, the substitution
    u = r + t^2, y = sign(t) * sqrt(f(u) / 7) traces the whole (single)
    component smoothly through its vertical tangent at (r + 1, 0).

    If height(x, y) is given, each point is lifted to z = height(x, y), giving
    the graph of that function over the curve with the same parametrisation.
    """
    f = lambda u: u**3 + 3 * u**2 + 6
    r = next(z.real for z in np.roots([1, 3, 0, 6]) if abs(z.imag) < 1e-9)
    t_max = np.sqrt(x_max - 1 - r)

    def point(t):
        u = r + t**2
        x, y = u + 1, np.sign(t) * np.sqrt(max(f(u), 0) / 7)
        return np.array([x, y, height(x, y) if height else 0])

    return ParametricFunction(point, t_range=[-t_max, t_max, t_step], **kwargs)


def split_curve(curve, pieces=200):
    """The curve cut into equal slices of its parametrisation, in reverse order.

    Splitting lets each slice take its own colour, and curves split the same
    way transform slice-to-slice. Reversed so the lower branch (traced first)
    is drawn on top of the upper branch. Each slice overlaps the next by half a
    slice to hide the joins. (Not round caps: Manim's Cairo camera never resets
    the line cap, so it would leak into every later stroke, e.g. dotting the
    outlines that Write traces.)
    """
    group = VGroup(*(
        VMobject().pointwise_become_partial(
            curve, i / pieces, min((i + 1.5) / pieces, 1)
        )
        for i in range(pieces)
    ))
    group.submobjects.reverse()
    return group


def ball(point, color, radius=0.08):
    """A flat dot shaded to look like a ball, drawn above flat mobjects.

    Dot3D won't do: the 3D camera always draws spheres beneath flat mobjects
    like curves, except while they animate, so they flicker in front. Register
    this with add_fixed_orientation_mobjects so it faces the camera.
    """
    return Dot(point, radius=radius, color=color).set_sheen(0.3, UL).set_z_index(1)


def pi_creature(mode, color, height):
    """A 3Blue1Brown pi creature (SVGs from CaftBotti/manim_pi_creatures),
    with its body, the largest path, recoloured."""
    creature = SVGMobject(ASSETS / "PiCreature" / f"PiCreatures_{mode}.svg", height=height)
    body = max(creature.submobjects, key=lambda path: path.width * path.height)
    body.set_fill(color, opacity=1)
    return creature


def curve_x(y):
    """The x-coordinate of the point at height y on the elliptic curve above."""
    roots = np.roots([1, 3, 0, 6 - 7 * y**2])
    return next(z.real for z in roots if abs(z.imag) < 1e-9) + 1


def curve_y(x):
    """The y-coordinate of the point at x on the upper branch of the curve."""
    u = x - 1
    return np.sqrt(max(u**3 + 3 * u**2 + 6, 0) / 7)


def colored_math_tex(*parts, **kwargs):
    """MathTex from parts, each a plain string (white) or a (tex, color) pair."""
    parts = [(p, WHITE) if isinstance(p, str) else p for p in parts]
    tex = MathTex(*(t for t, _ in parts), **kwargs)
    for part, (_, color) in zip(tex, parts):
        part.set_color(color)
    return tex


def colored_tex(*parts, **kwargs):
    """Tex from parts, each a plain string (white) or a (tex, color) pair."""
    parts = [(p, WHITE) if isinstance(p, str) else p for p in parts]
    tex = Tex(*(t for t, _ in parts), **kwargs)
    for part, (_, color) in zip(tex, parts):
        part.set_color(color)
    return tex


class Preliminaries(ThreeDSlide):
    def construct(self):
        self.wait_time_between_slides = 0.2
        self.next_slide()

        curve = elliptic_curve(color=BLUE, stroke_width=4)
        # Same curve continued off-screen, swapped in for a later zoom-out. Split
        # like the function graphs so it transforms into them slice-to-slice.
        long_curve = split_curve(elliptic_curve(x_max=12)).set_stroke(BLUE, width=4)

        title = Tex(r"\textbf{Elliptic Curve}", r" \textbf{Codes:}", font_size=48)
        title.to_corner(UL, buff=0.5)

        subtitle = Tex(
            "The Genus 1 Advantage", font_size=40, tex_to_color_map={"Genus 1": BLUE}
        )
        subtitle.to_edge(RIGHT, buff=0.3)

        speakers = Tex(
            r"Jamie Low Jiaqi \\ Oliver Ang Wei Jie",
            font_size=36,
            tex_environment="flushleft",
        )
        speakers.to_corner(DL, buff=0.5)

        logo = ImageMobject(ASSETS / "ntu_spms_logo_monochrome_inverted.png")
        logo.scale_to_fit_height(0.5).next_to(speakers, RIGHT, buff=0.4)

        self.play(Create(curve), run_time=3)
        self.remove(curve)
        self.add(long_curve)
        self.play(Write(title))
        self.play(Write(subtitle))
        self.play(Write(speakers))
        self.play(FadeIn(logo))
        self.next_slide()

        self.play(Unwrite(title[1]), Unwrite(subtitle), Unwrite(speakers), FadeOut(logo))
        self.next_slide()

        # Powers of y in yellow, powers of x in green, everything else white.
        # Split over two lines so it fits inside the curve.
        lhs = colored_math_tex(
            ("y^2", YELLOW), ("+", WHITE), ("a_1", WHITE), ("x", GREEN),
            ("y", YELLOW), ("+", WHITE), ("a_3", WHITE), ("y", YELLOW),
            font_size=36,
        )
        rhs = colored_math_tex(
            ("=", WHITE), ("x^3", GREEN), ("+", WHITE), ("a_2", WHITE),
            ("x^2", GREEN), ("+", WHITE), ("a_4", WHITE), ("x", GREEN),
            ("+", WHITE), ("a_6", WHITE),
            font_size=36,
        )
        rhs.next_to(lhs, DOWN, buff=0.25, aligned_edge=LEFT).shift(0.5 * RIGHT)
        # Centre in the region enclosed by the curve and the right edge.
        weierstrass = VGroup(lhs, rhs).move_to(
            [(curve.get_left()[0] + config.frame_width / 2) / 2, 0, 0]
        )

        self.play(Write(weierstrass))
        self.next_slide()

        field = Tex(r"over $\mathbb{F}_q$", font_size=32)
        field.next_to(title[0], DOWN, buff=0.2, aligned_edge=LEFT)
        self.play(Write(field))

        # Highlight the coefficients a_1, a_3, a_2, a_4, a_6.
        coefficients = [lhs[2], lhs[6], rhs[3], rhs[6], rhs[9]]
        for _ in range(2):
            self.play(*(Indicate(a, color=PINK) for a in coefficients))
        self.next_slide()

        # Highlight the leading terms y^2 and x^3.
        for _ in range(2):
            self.play(Indicate(lhs[0], color=PINK), Indicate(rhs[1], color=PINK))
        self.next_slide()

        # Zoomed out, the frame spans x in [-12, 12] and y in [-9, 9]. s converts
        # on-screen sizes to scene units, so labels look like unzoomed text.
        zoom_out = 1 / 3
        s = 1 / zoom_out

        # Point at infinity: a label right of each branch near the frame edge,
        # with an arrow from its top (bottom) parallel to the curve there.
        label_y = 6.5
        u = curve_x(label_y) - 1
        slope = (3 * u**2 + 6 * u) / (14 * label_y)  # dy/dx on the upper branch
        infinities = VGroup()
        for sign in (1, -1):
            infinity = MathTex(r"\infty", font_size=36 * s)
            infinity.next_to([curve_x(label_y), sign * label_y, 0], RIGHT, buff=0.3 * s)
            edge = infinity.get_top() if sign > 0 else infinity.get_bottom()
            start = edge + sign * 0.15 * s * UP
            arrow = Arrow(
                start,
                start + 0.42 * s * normalize(RIGHT + sign * slope * UP),
                buff=0,
                stroke_width=3,
                tip_length=0.12 * s,
                max_tip_length_to_length_ratio=0.5,
            )
            infinities.add(infinity, arrow)
        self.add(infinities)

        self.move_camera(zoom=zoom_out)
        self.next_slide()

        self.move_camera(zoom=1)

        genus = Tex("Projective curve of genus 1", font_size=36)
        genus.to_corner(DL, buff=0.5)

        self.play(Write(genus))
        self.next_slide()

        # Fixed in frame so it stays put once the camera starts rotating. Register
        # it with the camera, then remove it so the transform can bring it in.
        ag_title = Tex(r"\textbf{AG Code}", font_size=40).to_corner(UL, buff=0.3)
        self.add_fixed_in_frame_mobjects(ag_title)
        self.remove(ag_title)

        # Axes in the curve's own coordinates: one unit per scene unit, origin at
        # the scene origin, extending past the frame so they still look long
        # enough once the camera rotates. A negative z-index keeps them behind the
        # curve, even while they are being drawn. The y range is asymmetric, so
        # shift the axes' origin back onto the scene origin.
        axes = ThreeDAxes(
            x_range=[-8, 8, 1], y_range=[-2, 3, 1], z_range=[-3, 3, 1],
            x_length=16, y_length=5, z_length=6,
            axis_config={"tip_width": 0.2, "tip_height": 0.2},
        ).set_z_index(-1)
        axes.shift(-axes.c2p(0, 0, 0))

        self.remove(infinities)  # off-screen since the zoom-in
        self.play(
            Unwrite(weierstrass),
            Unwrite(field),
            Unwrite(genus),
            ReplacementTransform(title[0], ag_title),
        )
        self.play(Create(axes))
        self.next_slide()

        # Tilt into a 3D view; at this angle the curve stays clear of "AG Code".
        self.move_camera(phi=65 * DEGREES, theta=-60 * DEGREES, run_time=2)

        # Eight evaluation points at random places on the visible part of the
        # curve, one in each eighth of the curve's parameter range so they don't
        # bunch up (fixed seed, so every render picks the same points). The
        # parameter is t = sign(y) * sqrt(x - x_left), as in elliptic_curve.
        rng = np.random.default_rng(4310)
        x_left = curve.get_left()[0]
        t_max = np.sqrt(config.frame_width / 2 - x_left)
        ts = (np.arange(8) + rng.uniform(0.2, 0.8, 8)) / 8 * 2 * t_max - t_max
        points = VGroup(*(
            ball([x_left + t**2, np.sign(t) * curve_y(x_left + t**2), 0], YELLOW)
            for t in ts
        ))
        self.add_fixed_orientation_mobjects(*points)
        self.remove(*points)
        self.play(LaggedStart(*(GrowFromCenter(p) for p in points), lag_ratio=0.2))

        # Labels stacked in the bottom-left corner: the pole group (written later)
        # at the bottom, the evaluation group directly above it.
        label_buff = 0.15
        pole_divisor = VGroup(
            Tex("Pole budget: $k$", font_size=28),
            MathTex(r"G = n_1(Q_1) + \cdots + n_r(Q_r)", font_size=28),
        ).set_color(RED).arrange(DOWN, buff=label_buff, aligned_edge=LEFT)
        # Left of the default corner spot to clear the z-axis, and a little lower.
        pole_divisor.to_corner(DL, buff=0.5).shift(0.2 * LEFT + 0.15 * DOWN)

        divisor = VGroup(
            Tex("Evaluation points:", font_size=28),
            MathTex(r"D = (P_1) + \cdots + (P_n)", font_size=28),
        ).set_color(YELLOW).arrange(DOWN, buff=label_buff, aligned_edge=LEFT)
        divisor.next_to(pole_divisor, UP, buff=label_buff, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(divisor)
        self.remove(divisor)
        self.play(Write(divisor))
        self.next_slide()

        # Two pole points, both on the curve.
        poles = VGroup(*(
            ball(point, RED) for point in ([-2, curve_y(-2), 0], [0.73, 0.94105, 0])
        ))
        self.add_fixed_orientation_mobjects(*poles)
        self.remove(*poles)
        self.play(LaggedStart(*(GrowFromCenter(p) for p in poles), lag_ratio=0.2))

        self.add_fixed_in_frame_mobjects(pole_divisor)
        self.remove(pole_divisor)
        self.play(Write(pole_divisor))
        self.next_slide()
        # A function on the curve, drawn as its graph z = f(x, y) above the curve:
        # f always lies in [0.2, 0.9], so it is positive and bounded by 1.
        def f(x, y):
            return 0.55 + 0.35 * np.tanh(y) * np.cos(1.5 * x)

        # The graph of g over the curve, shaded light green (lower branch) to dark
        # green (upper branch) so crossing branches are easy to tell apart. One
        # path can only take one colour, so it is cut into short pieces. Fine
        # sampling keeps the spikes at the poles added later smooth.
        # The gradient spans only the on-screen part of the curve (x <= 4.2, the
        # same as the blue curve), so both ends of it are visible; beyond that
        # the colour stays at the end values.
        def graph_of(g, pieces=200, light="#E6F59A", dark="#216E53"):
            graph = elliptic_curve(x_max=12, height=g, t_step=0.002, stroke_width=4)
            # Pieces are uniform in the parameter t, which runs over
            # [-t_full, t_full]; the on-screen part is [-t_visible, t_visible].
            t_full = np.sqrt(12 - x_left)
            t_visible = np.sqrt(config.frame_width / 2 + 0.2 - x_left)
            group = split_curve(graph, pieces)
            for j, piece in enumerate(group):
                i = pieces - 1 - j  # split_curve reverses the slice order
                t = ((i + 0.5) / pieces * 2 - 1) * t_full
                alpha = np.clip((t + t_visible) / (2 * t_visible), 0, 1)
                color = interpolate_color(ManimColor(light), ManimColor(dark), alpha)
                piece.set_stroke(color, width=4)
            return group

        f_graph = graph_of(f)
        self.play(ReplacementTransform(long_curve.copy(), f_graph), run_time=2)

        # Dotted lines from each point on the curve up to the graph of func.
        def lines_to_graph(dots, color, func=f):
            return VGroup(*(
                DashedLine(
                    (p := dot.get_center()),
                    [p[0], p[1], func(p[0], p[1])],
                    dash_length=0.04,
                    dashed_ratio=0.5,
                    stroke_width=2,
                    color=color,
                )
                for dot in dots
            ))

        evaluations = lines_to_graph(points, YELLOW)
        pole_lines = lines_to_graph(poles, RED)
        self.play(LaggedStart(*(Create(line) for line in evaluations + pole_lines), lag_ratio=0.1))
        self.next_slide()

        # Give a base function poles at red points: for each (Q, c), add a term
        # that behaves like c / |P - Q|^2 near Q (spiking up for c > 0, down for
        # c < 0) but is negligible beyond |P - Q| ~ 1. Values are clipped to
        # [-6, 6], well beyond the frame.
        def with_poles(base, *poles_and_coefficients):
            def g(x, y):
                value = base(x, y)
                for q, c in poles_and_coefficients:
                    d2 = (x - q[0]) ** 2 + (y - q[1]) ** 2
                    value += c * np.exp(-4 * d2) / max(d2, 1e-9)
                return np.clip(value, -6, 6)
            return g

        q1, q2 = poles[0].get_center(), poles[1].get_center()

        # Animate f (graph and dotted lines) into a new function g.
        def morph_f_into(g):
            self.play(
                Transform(f_graph, graph_of(g)),
                Transform(evaluations, lines_to_graph(points, YELLOW, g)),
                Transform(pole_lines, lines_to_graph(poles, RED, g)),
                run_time=2,
            )

        # A pole at the right-hand red point...
        morph_f_into(with_poles(f, (q2, 0.02)))
        self.next_slide()

        # ...and a second pole at the left-hand red point.
        morph_f_into(with_poles(f, (q2, 0.02), (q1, 0.02)))
        self.next_slide()

        # A different function, taking both signs, with a pole at the
        # right-hand red point only.
        def g(x, y):
            return 0.7 * np.cos(1.3 * x + 1.2 * y)

        morph_f_into(with_poles(g, (q2, 0.02)))
        self.next_slide()

        # Another function with poles at the same points: a wave that changes
        # sign along the curve, with the pole at the right-hand point flipped.
        def h(x, y):
            return 0.8 * np.sin(1.8 * x - 0.6 * y)

        h_with_poles = with_poles(h, (q2, -0.02), (q1, 0.02))
        morph_f_into(h_with_poles)
        self.next_slide()

        # Top-right labels: the function space, with the codewords (written on
        # the next slide) below it. Two codeword lines (one Tex) keep it clear of
        # the z-axis. Placed together so the wider one sets the corner gap.
        function_space = Tex(
            r"Function space: $f \in \mathcal{L}(G)$", font_size=28, color=RED
        )
        codewords = Tex(
            r"Codewords: \\ $f(D) := (f(P_1), \hdots, f(P_n))$",
            font_size=28,
            color=YELLOW,
            tex_environment="flushleft",
        )
        VGroup(function_space, codewords).arrange(
            DOWN, buff=label_buff, aligned_edge=LEFT
        ).to_corner(UR, buff=0.3)
        self.add_fixed_in_frame_mobjects(function_space)
        self.remove(function_space)
        self.play(Write(function_space))
        self.next_slide()

        # Evaluate: each evaluation point rises (or falls) along its dotted line
        # to the graph.
        targets = [
            [x, y, h_with_poles(x, y)] for x, y, _ in (p.get_center() for p in points)
        ]
        self.play(
            *(point.animate.move_to(target) for point, target in zip(points, targets)),
            run_time=2,
        )

        self.add_fixed_in_frame_mobjects(codewords)
        self.remove(codewords)
        self.play(Write(codewords))

        code_name = colored_math_tex(
            (r"C_{\mathcal{L}}(", WHITE), ("D", YELLOW), (",", WHITE), ("G", RED),
            (")", WHITE),
            font_size=36,
        ).next_to(ag_title, DOWN, buff=0.2, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(code_name)
        self.remove(code_name)
        self.play(Write(code_name))
        self.next_slide()

        # n yellow (evaluation points), k red (pole budget), EC blue (the curve).
        parameters = colored_math_tex(
            ("[", WHITE), ("n", YELLOW), (",", WHITE), ("k", RED),
            (r",\geq\!", WHITE), ("n", YELLOW), ("-", WHITE), ("k", RED),
            ("]", WHITE), (r"\text{ for }", WHITE), (r"\text{EC}", BLUE),
            (r"\text{s}", WHITE),
            font_size=28,
        )
        parameters.next_to(codewords, DOWN, buff=label_buff, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(parameters)
        self.remove(parameters)
        self.play(Write(parameters))
        self.next_slide()

        # Clear everything but the blue curve, then return to the flat view.
        self.play(
            *(Unwrite(text) for text in (
                ag_title, code_name, divisor, pole_divisor, function_space,
                codewords, parameters,
            )),
            *(Uncreate(mob) for mob in (
                axes, f_graph, evaluations, pole_lines, *points, *poles,
            )),
            run_time=2,
        )
        self.move_camera(phi=0, theta=-90 * DEGREES, run_time=2)
        self.next_slide()

        group_law_title = Tex(r"\textbf{Elliptic Group Law}", font_size=40)
        group_law_title.to_corner(UL, buff=0.3)
        self.play(Write(group_law_title))
        self.next_slide()

        # Chord construction. The curve is 7y^2 = x^3 - 3x + 8 (no x^2 term), so
        # the three x-coordinates where the line y = mx + b meets it sum to
        # 7m^2, which gives the third point R. P (+) Q is R reflected in the
        # x-axis.
        p = np.array([-2.4, curve_y(-2.4), 0])
        q = np.array([0.6, -curve_y(0.6), 0])
        m = (q[1] - p[1]) / (q[0] - p[0])
        x_r = 7 * m**2 - p[0] - q[0]
        r = np.array([x_r, p[1] + m * (x_r - p[0]), 0])
        p_plus_q = np.array([r[0], -r[1], 0])

        def labelled_point(point, tex, direction, buff=0.15):
            dot = Dot(point, radius=0.08, color=WHITE).set_z_index(1)
            label = MathTex(tex, font_size=36).next_to(dot, direction, buff=buff)
            return dot, label

        p_dot, p_label = labelled_point(p, "P", UP)
        q_dot, q_label = labelled_point(q, "Q", DOWN)
        self.play(GrowFromCenter(p_dot), GrowFromCenter(q_dot))
        self.play(Write(p_label), Write(q_label))
        self.next_slide()

        # The chord through P and Q, running a little past P and R.
        direction = normalize(r - p)
        chord = Line(p - 0.5 * direction, r + 0.5 * direction, color=GREEN)
        r_dot, r_label = labelled_point(r, "R", DL, buff=0.03)
        self.play(Create(chord))
        self.play(GrowFromCenter(r_dot))
        self.play(Write(r_label))
        self.next_slide()

        # Reflect R in the x-axis: the vertical line from R meets the curve again
        # at P (+) Q.
        vertical = Line(r, p_plus_q, color=GREEN)
        sum_dot, sum_label = labelled_point(p_plus_q, r"P \oplus Q", LEFT)
        self.play(Create(vertical))
        self.play(GrowFromCenter(sum_dot))
        self.play(Write(sum_label))
        self.next_slide()

        # Clear the construction, then slide P and Q along their branches to
        # x = 1, one directly above the other.
        self.play(
            Uncreate(chord), Uncreate(vertical),
            FadeOut(r_dot), FadeOut(sum_dot), Unwrite(r_label), Unwrite(sum_label),
        )
        p_x, q_x = ValueTracker(p[0]), ValueTracker(q[0])
        p_dot.add_updater(lambda d: d.move_to([p_x.get_value(), curve_y(p_x.get_value()), 0]))
        q_dot.add_updater(lambda d: d.move_to([q_x.get_value(), -curve_y(q_x.get_value()), 0]))
        # Labels move off to the right as the points arrive, clear of the
        # vertical line drawn through them next.
        p_label.add_updater(lambda l: l.next_to(p_dot, UR, buff=0.05))
        q_label.add_updater(lambda l: l.next_to(q_dot, DR, buff=0.05))
        self.play(p_x.animate.set_value(1), q_x.animate.set_value(1), run_time=2)
        for mob in (p_dot, q_dot, p_label, q_label):
            mob.clear_updaters()
        self.next_slide()

        # The vertical line through P and Q, drawn upwards across the frame. Then
        # it is swapped for one reaching well past the frame even after zooming
        # out; the swap is invisible since they agree on screen.
        vertical_pq = Line([1, -3.2, 0], [1, 3.2, 0], color=GREEN)
        self.play(Create(vertical_pq), run_time=1.5)
        self.remove(vertical_pq)
        vertical_pq = Line([1, -40, 0], [1, 40, 0], color=GREEN)
        self.add(vertical_pq)
        self.next_slide()

        # Zoom out to 1/3 (frame y in [-9, 9]). The label sits at the top, where
        # the line was heading, left of the line. It is added beforehand off the
        # top of the current frame, so it comes into view as the camera pulls
        # back. Sized x3 so it reads at the usual size once zoomed.
        zoom_out = 1 / 3
        s = 1 / zoom_out
        infinity_label = MathTex(r"\infty = P \oplus Q", font_size=36 * s)
        infinity_label.next_to([1, 7.9, 0], LEFT, buff=0.15 * s)
        self.add(infinity_label)
        self.move_camera(zoom=zoom_out, run_time=2)
        self.next_slide()

        # The inverse of P, left of the line just below Q (and below the loop).
        inverse_label = MathTex(r"Q = \ominus P", font_size=36 * s)
        inverse_label.next_to([1, q_dot.get_y() - 1.3, 0], LEFT, buff=0.15 * s)
        self.play(Write(inverse_label))
        self.next_slide()

        # Back to the normal view, dropping the labels on the way. Swap the long
        # line for its on-screen part so Uncreate is visible throughout.
        self.move_camera(
            zoom=1,
            run_time=2,
            added_anims=[FadeOut(infinity_label), FadeOut(inverse_label)],
        )
        self.remove(vertical_pq)
        vertical_pq = Line([1, -3.2, 0], [1, 3.2, 0], color=GREEN)
        self.add(vertical_pq)
        self.play(Uncreate(vertical_pq))

        # Merge P and Q into a single point P on the upper branch.
        p = np.array([2.5, curve_y(2.5), 0])
        merged_label = MathTex("P", font_size=36).next_to(p, UL, buff=0.15)
        self.play(
            p_dot.animate.move_to(p),
            q_dot.animate.move_to(p),
            ReplacementTransform(VGroup(p_label, q_label), merged_label),
            run_time=2,
        )
        self.remove(q_dot)
        p_label = merged_label
        self.next_slide()

        # Tangent construction. Implicit differentiation of 7y^2 = x^3 - 3x + 8
        # gives the slope m = (3x^2 - 3) / (14y). The tangent meets the curve
        # twice at P, so x_R = 7m^2 - 2 x_P. [2]P is R reflected in the x-axis.
        m = (3 * p[0] ** 2 - 3) / (14 * p[1])
        x_r = 7 * m**2 - 2 * p[0]
        r = np.array([x_r, p[1] + m * (x_r - p[0]), 0])
        double_p = np.array([r[0], -r[1], 0])

        # The tangent, drawn from beyond P to beyond R.
        direction = normalize(r - p)
        tangent = Line(p - 1.6 * direction, r + 0.6 * direction, color=GREEN)
        r_dot, r_label = labelled_point(r, "R", DR, buff=0.05)
        r_label.shift(0.15 * LEFT)  # tuck it in under R
        self.play(Create(tangent))
        self.play(GrowFromCenter(r_dot))
        self.play(Write(r_label))

        # Reflect R in the x-axis: the vertical line from R meets the curve again
        # at [2]P.
        vertical = Line(r, double_p, color=GREEN)
        double_dot, double_label = labelled_point(double_p, "[2]P", UP)
        self.play(Create(vertical))
        self.play(GrowFromCenter(double_dot))
        self.play(Write(double_label))
        self.next_slide()

        # R is the inverse of [2]P.
        minus_double_label = MathTex("[-2]P", font_size=36)
        minus_double_label.next_to(r_dot, DR, buff=0.05).shift(0.15 * LEFT)
        self.play(ReplacementTransform(r_label, minus_double_label))
        self.next_slide()

        # Clear all but P, then slide P back along the upper branch to the
        # leftmost point of the curve (a Weierstrass point: y = 0, so P = -P).
        self.play(
            Uncreate(tangent), Uncreate(vertical),
            FadeOut(r_dot), FadeOut(double_dot),
            Unwrite(minus_double_label), Unwrite(double_label),
        )
        x_weierstrass = curve_x(0)
        p_x = ValueTracker(p[0])
        p_dot.add_updater(lambda d: d.move_to([p_x.get_value(), curve_y(p_x.get_value()), 0]))
        p_label.add_updater(lambda l: l.next_to(p_dot, UL, buff=0.15))
        self.play(p_x.animate.set_value(x_weierstrass), run_time=2)
        p_dot.clear_updaters()
        p_label.clear_updaters()
        self.next_slide()

        # The tangent at a Weierstrass point is vertical, so it meets the curve
        # again only at infinity: [2]P = infinity. Drawn upwards across the
        # frame while the title becomes that equation, top right of the line.
        vertical_tangent = Line(
            [x_weierstrass, -3.2, 0], [x_weierstrass, 3.2, 0], color=GREEN
        )
        doubling_infinity = MathTex(r"\infty = [2]P", font_size=36)
        doubling_infinity.next_to([x_weierstrass, 2.6, 0], RIGHT, buff=0.15)
        self.play(
            Create(vertical_tangent),
            ReplacementTransform(group_law_title, doubling_infinity),
            run_time=2,
        )
        self.next_slide()

        # Clear the stage.
        self.play(
            Unwrite(doubling_infinity), Unwrite(p_label),
            Uncreate(vertical_tangent), Uncreate(p_dot), Uncreate(long_curve),
            run_time=2,
        )
        # Unwrite's last frame leaves specks of the glyphs; the text is only
        # removed after it. A short wait renders one clean frame to end on.
        self.wait(0.1)


class Heart(Slide):
    def construct(self):
        self.wait_time_between_slides = 0.2

        # As before: D, the P_i and n yellow; G, the Q_i and k red. The new
        # divisor G' and its point Q are pink. Multiplicities stay white, as
        # does k inside the theorem statement.
        y, r, p = YELLOW, RED, PINK

        # Long statements are split into separately built lines: line breaks
        # inside a multi-part Tex don't land reliably.
        recall = VGroup(
            colored_tex(
                r"Recall that the code $C = C_{\mathcal{L}}(", ("D", y), ",",
                ("G", r), ")$",
                font_size=28,
            ),
            colored_tex(
                r"has parameters $[", ("n", y), ",", ("k", r), r",\geq\!", ("n", y),
                "-", ("k", r), "]$",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.15)
        divisor_d = colored_math_tex(
            ("D", y), "=(", ("P_1", y), r")+\cdots+(", ("P_n", y), ")",
            font_size=28,
        )
        divisor_g = colored_tex(
            "$", ("G", r), "=n_1(", ("Q_1", r), r")+\cdots+n_r(", ("Q_r", r),
            r")$ with $\sum_{i=1}^r n_i=", ("k", r), "$",
            font_size=28,
        )
        heading = Tex(r"\textbf{Two-Point Normal Form}", font_size=36)
        normal_form = VGroup(
            colored_tex(
                "Up to coordinate scaling, $", ("G", r), "$ may be replaced by",
                font_size=28,
            ),
            colored_tex(
                "$", ("G'", p), r"=(k-1)(\infty)+(", ("Q", p), ")$", ", where",
                font_size=28,
            ),
            colored_tex(
                "$", ("Q", p), "=[n_1]", ("Q_1", r), r"\oplus\cdots\oplus[n_r]",
                ("Q_r", r), "$.",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.15)

        recall_block = VGroup(recall, divisor_d, divisor_g).arrange(DOWN, buff=0.25)
        VGroup(recall_block, heading, normal_form).arrange(DOWN, buff=0.6)

        for tex in recall_block:
            self.play(Write(tex))
        self.next_slide()

        self.play(Write(heading))
        self.play(Write(normal_form))
        self.next_slide()

        # G' takes the place of G in the recall block...
        g_prime_line = normal_form[1]
        g_prime = VGroup(*g_prime_line[:-1])  # without ", where"
        where = g_prime_line[-1]
        self.play(
            FadeOut(divisor_g), FadeOut(where), g_prime.animate.move_to(divisor_g)
        )

        # ...then the theorem goes, and G' becomes the new G.
        new_g = colored_tex(
            "$", ("G", r), "=(", ("k", r), r"-1)(\infty)+(", ("Q", r), ")$",
            font_size=28,
        ).move_to(g_prime)
        self.play(
            Unwrite(heading),
            Unwrite(normal_form[0]), Unwrite(normal_form[2]),
            ReplacementTransform(g_prime, new_g),
        )
        self.next_slide()

        # The criterion, in the space the normal form left behind.
        criterion_heading = Tex(
            r"\textbf{Subset-Sum Criterion} \\ \textbf{for MDS Property}", font_size=36
        )
        criterion = VGroup(
            colored_tex("$C$ is MDS if and only if no $", ("k", r), "$ points", font_size=28),
            colored_tex(
                "$", ("P_{i_1}", y), r",\hdots,", ("P_{i_k}", y), "$ sum to $",
                ("Q", r), "$",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.15)
        VGroup(criterion_heading, criterion).arrange(DOWN, buff=0.4).next_to(
            recall_block, DOWN, buff=0.6
        )
        footnote = Tex(
            "This is a special case of the full theorem, which can be found in our report.",
            font_size=18,
        ).to_edge(DOWN, buff=0.3)

        self.play(Write(criterion_heading))
        self.play(Write(criterion))
        self.play(Write(footnote))
        self.next_slide()

        # Clear the recall block and footnote, then condense the criterion into
        # one line at the top, with the definition of iso-duality beneath it.
        self.play(
            Unwrite(recall), Unwrite(divisor_d), Unwrite(new_g), Unwrite(footnote)
        )
        summary_size = 22
        criterion_line = colored_tex(
            r"\textbf{Subset-Sum Criterion:} $C$ is MDS $\iff ", ("P_{i_1}", y),
            r"\oplus\cdots\oplus ", ("P_{i_k}", y), r"\neq ", ("Q", r), "$",
            font_size=summary_size,
        )
        iso_dual = Tex(
            r"$C$ is \textbf{iso-dual} if $C^\perp$ is a scaling of $C$",
            font_size=summary_size,
        )
        VGroup(criterion_line, iso_dual).arrange(DOWN, buff=0.15).to_edge(UP, buff=0.4)
        self.play(ReplacementTransform(VGroup(criterion_heading, criterion), criterion_line))
        self.play(Write(iso_dual))
        self.next_slide()

        # The construction, one clause per line so each fits the frame, in three
        # spaced groups: choosing G, choosing D, and the conclusion. k is white
        # here. Kept compact to leave room for the argument below it.
        size = 26
        iso_dual_heading = Tex(r"\textbf{Iso-Dual MDS Elliptic Codes}", font_size=32)
        construction = VGroup(
            VGroup(
                colored_tex(
                    "Choose a point $", ("Q", r), "$ of order 2 and $k$ even.",
                    font_size=size,
                ),
                colored_tex(
                    "Set $", ("G", r), r"=(k-1)(\infty)+(", ("Q", r), ")$.",
                    font_size=size,
                ),
            ).arrange(DOWN, buff=0.12),
            VGroup(
                colored_tex(
                    "Choose points $", ("P_1", y), r",\hdots,", ("P_k", y),
                    "$ of odd order.*",
                    font_size=size,
                ),
                colored_tex(
                    "Set $", ("D", y), r"=\sum_{i=1}^{k}\bigl((", ("Q", r),
                    r"\oplus ", ("P_i", y), ")+(", ("Q", r), r"\ominus ",
                    ("P_i", y), r")\bigr)$.",
                    font_size=size,
                ),
            ).arrange(DOWN, buff=0.12),
            colored_tex(
                r"Then $C_{\mathcal{L}}(", ("D", y), ",", ("G", r),
                ")$ is a $[2k,k,k+1]$ ", "iso-dual", " MDS", " code.",
                font_size=size,
            ),
        ).arrange(DOWN, buff=0.3)
        VGroup(iso_dual_heading, construction).arrange(DOWN, buff=0.3).next_to(
            iso_dual, DOWN, buff=0.35
        )
        conditions = Tex(
            r"*The points must be distinct, not $\infty$, and no two are inverses.",
            font_size=18,
        ).to_edge(DOWN, buff=0.3)
        self.play(Write(iso_dual_heading))
        self.play(Write(construction[0]))
        self.next_slide()
        self.play(Write(construction[1]))
        self.next_slide()
        self.play(Write(construction[2]))
        self.play(Write(conditions))
        self.next_slide()

        # Why it is MDS: the subset-sum criterion.
        self.play(Indicate(construction[2][-2], color=PINK), Unwrite(conditions))  # "MDS"
        self.next_slide()

        self.play(Indicate(criterion_line, color=PINK, scale_factor=1.05))
        self.next_slide()

        # Summing k of the points Q + P_i: Q has order 2 and k is even, so
        # [k]Q is infinity, and the odd-order P_i can't sum to Q. Built as one
        # MathTex for consistent spacing, then revealed segment by segment, each
        # derived from a copy of the previous one.
        chain = colored_math_tex(
            r"\bigoplus_{i=1}^k (", ("Q", r), r"\oplus ", ("P_i", y), ")",  # 0-4
            "=",  # 5
            "[k]", ("Q", r), r"\oplus\bigoplus_{i=1}^k ", ("P_i", y),  # 6-9
            "=",  # 10
            r"\bigoplus_{i=1}^k ", ("P_i", y),  # 11-12
            r"\neq",  # 13
            ("Q", r),  # 14
            font_size=28,
        ).to_edge(DOWN, buff=0.5)
        segments = [chain[0:5], chain[6:10], chain[11:13], chain[14:15]]
        relations = [chain[5], chain[10], chain[13]]

        self.play(Write(segments[0]))
        for relation, previous, segment in zip(relations, segments, segments[1:]):
            self.next_slide()
            self.play(Write(relation))
            self.play(TransformFromCopy(previous, segment) if segment != segments[-1] else Write(segment))
        self.next_slide()

        # So no k points sum to Q: the code is MDS.
        self.play(Indicate(construction[2][-2], color=PINK))  # "MDS"
        self.next_slide()

        # Now iso-duality: the plan of attack.
        self.play(Unwrite(chain))
        self.play(Indicate(construction[2][-3], color=PINK))  # "iso-dual"
        idea = VGroup(
            colored_tex(
                r"\textit{Idea:} Show that $C_{\mathcal{L}}(", ("D", y), ",", ("G", r),
                r")^\perp=C_{\mathcal{L}}(", ("D", y), ",", ("H", p), ")$",
                font_size=28,
            ),
            colored_tex(
                "where the points in $", ("H", p), "$ sum to $", ("Q", r), "$.",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.15).to_edge(DOWN, buff=0.6)
        self.play(Write(idea))
        self.next_slide()

        self.play(*(
            Unwrite(mob)
            for mob in (criterion_line, iso_dual, iso_dual_heading, construction, idea)
        ))
        # Unwrite's last frame leaves specks; end on a clean frame.
        self.wait(0.1)


class Isogenies(Slide):
    def construct(self):
        self.wait_time_between_slides = 0.2

        title = Tex(r"\textbf{Isogenies}", font_size=40).to_edge(UP, buff=0.4)

        # E_2 is the curve from the start, 7y^2 = x^3 - 3x + 8. E_1 is the image
        # of a degree-4 isogeny: two 2-isogenies, each by Velu's formula
        #     y^2 = x(x^2 + ax + b)  ->  y^2 = x((x - a)^2 - 4b),
        # whose kernel is the 2-torsion point moved to (0, 0). The first uses
        # E_2's only real 2-torsion point and gives a curve with two components
        # whose oval and branch almost touch; the second uses that curve's
        # rightmost root, which keeps two components and spreads them apart.
        r = curve_x(0)  # E_2's real root; in x - r, E_2 is x(x^2 + ax + b)
        a, b = 3 * r, 3 * r**2 - 3
        low, high = a - 2 * np.sqrt(b), a + 2 * np.sqrt(b)  # middle curve's roots
        a2, b2 = 2 * high - low, (high - low) * high  # moving `high` to 0
        c1, c2 = a2 - 2 * np.sqrt(b2), a2 + 2 * np.sqrt(b2)  # E_1: x(x-c1)(x-c2)

        # Each graph on faint axes at equal x and y scales, so shapes are true.
        def curve_graph(f, x_range, y_range, height=2.8):
            unit = height / (y_range[1] - y_range[0])
            axes = Axes(
                x_range=[*x_range, 1], y_range=[*y_range, 1],
                x_length=unit * (x_range[1] - x_range[0]), y_length=height,
                tips=False, axis_config={"stroke_opacity": 0.4, "include_ticks": False},
            )
            graph = axes.plot_implicit_curve(f, color=BLUE, stroke_width=3)
            return VGroup(axes, graph)

        e2 = curve_graph(lambda x, y: 7 * y**2 - (x**3 - 3 * x + 8), (-3, 4), (-3, 3))
        e1 = curve_graph(
            lambda x, y: 7 * y**2 - x * (x - c1) * (x - c2), (-3, 25), (-14, 14)
        )
        VGroup(e1, e2).arrange(RIGHT, buff=1.0).move_to(0.2 * UP)
        e1_label = MathTex(r"\mathcal{E}_1", font_size=36).next_to(e1, DOWN, buff=0.2)
        e2_label = MathTex(r"\mathcal{E}_2", font_size=36).next_to(e2, DOWN, buff=0.2)
        arrow = Arrow(e1.get_right(), e2.get_left(), buff=0.15)
        arrow_label = MathTex(r"\varphi", font_size=36).next_to(arrow, UP, buff=0.1)

        self.play(Write(title))
        self.play(Create(e1), Create(e2))
        self.play(Write(e1_label), Write(e2_label))
        self.play(GrowArrow(arrow), Write(arrow_label))

        description = Tex(
            "Curve morphisms that preserve the group structure", font_size=28
        ).to_edge(DOWN, buff=0.4)
        self.play(Write(description))
        self.next_slide()

        # The arrow E_1 -> E_2 is the dual of that isogeny. (The isogeny itself
        # wraps E_2's real points 4 times around E_1's branch, so it can't send
        # two points and their sum into E_1's window; the dual is one-to-one on
        # the branch.) Each step is the dual of a Velu 2-isogeny, with 7y^2:
        #     7y^2 = x((x - a)^2 - 4b)  ->  7y^2 = x(x^2 + ax + b),
        #     (x, y) -> (7y^2 / 4x^2, y(a^2 - 4b - x^2) / 8x^2).
        def velu_dual(a, b, x, y):
            return 7 * y**2 / (4 * x**2), y * (a**2 - 4 * b - x**2) / (8 * x**2)

        def isogeny_dual(x, y):
            x, y = velu_dual(a2, b2, x, y)  # to the middle curve, `high` at 0
            x, y = velu_dual(a, b, x + high, y)  # to E_2, its real root at 0
            return x + r, y

        # Chord construction on 7y^2 = x^3 + Ax^2 + ...: the line y = mx + c
        # meets the curve where the x-coordinates sum to 7m^2 - A.
        def chord_third_point(p, q, x2_coefficient):
            m = (q[1] - p[1]) / (q[0] - p[0])
            x = 7 * m**2 - x2_coefficient - p[0] - q[0]
            return np.array([x, p[1] + m * (x - p[0])])

        def e1_point(x, sign):
            return np.array([x, sign * np.sqrt(x * (x - c1) * (x - c2) / 7)])

        # P on the oval, Q on the branch; R and P (+) Q land on the oval.
        p, q = e1_point(9.8, 1), e1_point(21.8, -1)
        r1 = chord_third_point(p, q, -(c1 + c2))
        sum1 = r1 * [1, -1]
        # Their images on E_2, and the chord construction there.
        p2, q2, sum2 = (np.array(isogeny_dual(*point)) for point in (p, q, sum1))
        r2 = chord_third_point(p2, q2, 0)

        def dot(graph, point):
            return Dot(graph[0].c2p(*point), radius=0.06, color=WHITE).set_z_index(1)

        def construction(graph, p, q, r, total):
            """The chord through p and q (a little past both ends) and the
            vertical line from the third point r to the sum."""
            ends = sorted((p, q, r), key=lambda point: point[0])
            direction = normalize(ends[-1] - ends[0])
            span = np.linalg.norm(ends[-1] - ends[0])
            chord = Line(
                graph[0].c2p(*(ends[0] - 0.06 * span * direction)),
                graph[0].c2p(*(ends[-1] + 0.06 * span * direction)),
                color=GREEN, stroke_width=3,
            )
            vertical = Line(graph[0].c2p(*r), graph[0].c2p(*total), color=GREEN, stroke_width=3)
            return chord, vertical

        p_dot, q_dot, r1_dot, sum1_dot = (dot(e1, point) for point in (p, q, r1, sum1))
        chord1, vertical1 = construction(e1, p, q, r1, sum1)
        self.play(GrowFromCenter(p_dot), GrowFromCenter(q_dot))
        self.play(Create(chord1))
        self.play(GrowFromCenter(r1_dot))
        self.play(Create(vertical1))
        self.play(GrowFromCenter(sum1_dot))
        self.next_slide()

        # Map P, Q and P (+) Q across, then repeat the construction on E_2: it
        # lands on the image of P (+) Q.
        moving = [d.copy() for d in (p_dot, q_dot, sum1_dot)]
        self.play(*(
            d.animate.move_to(e2[0].c2p(*target))
            for d, target in zip(moving, (p2, q2, sum2))
        ), run_time=2)
        r2_dot = dot(e2, r2)
        chord2, vertical2 = construction(e2, p2, q2, r2, sum2)
        self.play(Create(chord2))
        self.play(GrowFromCenter(r2_dot))
        self.play(Create(vertical2))
        self.next_slide()

        # Clear both constructions. Then mark E_1's Weierstrass points (its
        # 2-torsion points, y = 0) as the map becomes doubling, [2], from E_1 to
        # itself: E_2 is replaced by a copy of E_1 in the same place.
        self.play(*(
            Uncreate(mob)
            for mob in (
                p_dot, q_dot, r1_dot, sum1_dot, *moving, r2_dot,
                chord1, vertical1, chord2, vertical2,
            )
        ))
        weierstrass = [dot(e1, (x, 0)) for x in (0, c1, c2)]
        # The point at infinity completes the kernel of [2]; like its image
        # later, it sits above the graph, flush with the right edge.
        def infinity_point(graph):
            point = Dot(radius=0.06, color=WHITE).set_z_index(1)
            point.next_to(graph, UP, buff=0.25).align_to(graph, RIGHT)
            label = MathTex(r"\infty", font_size=36).next_to(point, LEFT, buff=0.15)
            return point, label

        e1_infinity, e1_infinity_label = infinity_point(e1)
        kernel = [*weierstrass, e1_infinity]
        doubling = MathTex("[2]", font_size=36).move_to(arrow_label)
        target = e1.copy().align_to(e2, LEFT).set_y(e2.get_y())
        target_label = MathTex(r"\mathcal{E}_1", font_size=36).next_to(target, DOWN, buff=0.2)
        self.play(
            *(GrowFromCenter(d) for d in kernel),
            Write(e1_infinity_label),
            ReplacementTransform(arrow_label, doubling),
            ReplacementTransform(e2, target),
            ReplacementTransform(e2_label, target_label),
        )
        self.next_slide()

        # [2] sends the whole kernel to the identity: copies of all four points
        # merge into the point at infinity of the target copy of E_1.
        infinity, infinity_label = infinity_point(target)
        images = [d.copy() for d in kernel]
        self.play(*(d.animate.move_to(infinity) for d in images), run_time=2)
        self.remove(*images)
        self.add(infinity)
        self.play(Write(infinity_label))
        self.next_slide()

        # Going back: the arrow turns round and [2] becomes [2]^{-1}.
        inverse = MathTex("[2]^{-1}", font_size=36).next_to(arrow, UP, buff=0.1)
        self.play(Rotate(arrow, PI), ReplacementTransform(doubling, inverse))
        self.next_slide()

        # The captions and description give way to the point of it all.
        action = VGroup(
            colored_tex(
                r"Isogenies can act on $C_{\mathcal{L}}(", ("D", YELLOW), ",",
                ("G", RED), ")$ by pushing forward ($", (r"\varphi_*", BLUE), "$)",
                font_size=28,
            ),
            colored_tex(
                "or pulling back ($", (r"\varphi^*", BLUE), "$) points and functions",
                font_size=28,
            ),
        ).arrange(DOWN, buff=0.15).next_to(VGroup(e1, target), DOWN, buff=0.5).set_x(0)
        self.play(ReplacementTransform(VGroup(e1_label, target_label, description), action))
        self.next_slide()

        # Clear the graphs (with their points) and the arrow.
        self.play(*(
            FadeOut(mob)
            for mob in (
                e1, target, *kernel, e1_infinity_label, infinity, infinity_label,
                arrow, inverse,
            )
        ))

        # Commutative diagram, between the title and the sentence below:
        #   E_2 --f--> F_q
        #    ^        /
        #  phi   f o phi
        #    |    /
        #   E_1
        e2_node = MathTex(r"\mathcal{E}_2", font_size=40).move_to([-1.5, 1.3, 0])
        field_node = MathTex(r"\mathbb{F}_q", font_size=40).move_to([1.5, 1.3, 0])
        e1_node = MathTex(r"\mathcal{E}_1", font_size=40).move_to([-1.5, -1.0, 0])

        def diagram_arrow(start, end, label_parts, side):
            arrow = Arrow(
                start.get_center(), end.get_center(), buff=0.4,
                stroke_width=3, tip_length=0.18, max_tip_length_to_length_ratio=0.2,
            )
            label = colored_math_tex(*label_parts, font_size=32)
            label.next_to(arrow.get_center(), side, buff=0.15)
            return arrow, label

        phi = (r"\varphi", BLUE)
        f_arrow, f_label = diagram_arrow(e2_node, field_node, ["f"], UP)
        phi_arrow, phi_label = diagram_arrow(e1_node, e2_node, [phi], LEFT)
        pullback_arrow, pullback_label = diagram_arrow(
            e1_node, field_node, [r"f\circ", phi], DR
        )

        self.play(Write(e2_node), Write(field_node))
        self.play(GrowArrow(f_arrow), Write(f_label))
        self.next_slide()

        self.play(Write(e1_node))
        self.play(GrowArrow(phi_arrow), Write(phi_label))
        self.next_slide()

        self.play(GrowArrow(pullback_arrow), Write(pullback_label))
        self.next_slide()

        # Clear the diagram and title; the sentence moves up to head the slide.
        diagram = VGroup(
            e2_node, field_node, e1_node, f_arrow, f_label, phi_arrow, phi_label,
            pullback_arrow, pullback_label,
        )
        self.play(FadeOut(diagram), Unwrite(title))
        self.play(action.animate.to_edge(UP, buff=0.4))

        # The original code, then the four codes of Section 7 of the report:
        # pullbacks (phi^*) with D and G on E_2, pushforwards (phi_*) with D and
        # G on E_1. The isogenies are blue.
        y, r = YELLOW, RED
        up, down = (r"\varphi^*", BLUE), (r"\varphi_*", BLUE)
        original = colored_math_tex(
            r"C_{\mathcal{L}}(", ("D", y), ",", ("G", r), r")=\{f(", ("D", y),
            r"):f\in\mathcal{L}(", ("G", r), r")\}",
            font_size=28,
        )
        codes = [
            (r"C_{\mathcal{L}}(", up, ("D", y), ",", up, ("G", r), r")=\{f(", up,
             ("D", y), r"):f\in\mathcal{L}(", up, ("G", r), r")\}"),
            (r"C_{\mathcal{L}}(", down, ("D", y), ",", down, ("G", r), r")=\{f(", down,
             ("D", y), r"):f\in\mathcal{L}(", down, ("G", r), r")\}"),
            (up, r"C_{\mathcal{L}}(", ("D", y), ",", ("G", r), r")=\{(", up, "f)(", up,
             ("D", y), r"):f\in\mathcal{L}(", ("G", r), r")\}"),
            (down, r"C_{\mathcal{L}}(", ("D", y), ",", ("G", r), r")=\{(", down, "f)(",
             down, ("D", y), r"):f\in\mathcal{L}(", ("G", r), r")\}"),
        ]
        code_list = VGroup(*(
            VGroup(
                Tex(f"({i})", font_size=28),
                colored_math_tex(*parts, font_size=28),
            ).arrange(RIGHT, buff=0.25)
            for i, parts in enumerate(codes, start=1)
        )).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        VGroup(original, code_list).arrange(DOWN, buff=0.9).next_to(
            action, DOWN, buff=0.5
        ).set_x(0)
        transform_arrow = Arrow(
            original.get_bottom(), [0, code_list.get_top()[1], 0], buff=0.15,
            stroke_width=3, tip_length=0.12, max_tip_length_to_length_ratio=0.3,
        )
        self.play(Write(original))
        self.play(GrowArrow(transform_arrow))
        for item in code_list:
            self.play(TransformFromCopy(original, item))
        self.next_slide()

        # Drop the original code; the list moves up under the sentence, and the
        # experiment is introduced at the bottom.
        self.play(FadeOut(original), FadeOut(transform_arrow))
        self.play(code_list.animate.next_to(action, DOWN, buff=0.6).set_x(code_list.get_x()))
        experiment = VGroup(
            Tex(
                "We asked Claude to write SageMath code to analyse",
                font_size=26,
            ),
            colored_tex(
                "these four constructions for some $", (r"\varphi", BLUE), "$, $",
                ("D", YELLOW), "$, and $", ("G", RED), "$",
                font_size=26,
            ),
        ).arrange(DOWN, buff=0.12).to_edge(DOWN, buff=0.5)
        self.play(Write(experiment))
        self.next_slide()

        # Findings, one per construction, each written under its item.
        self.play(Unwrite(experiment))
        self.play(code_list[3].animate.shift(0.3 * DOWN))
        nonlinear = Tex("produces nonlinear codes", font_size=26).next_to(
            code_list[3][1], DOWN, buff=0.15, aligned_edge=LEFT
        )
        self.play(Write(nonlinear))
        self.next_slide()

        self.play(FadeOut(code_list[3]), FadeOut(nonlinear))
        self.play(code_list[2].animate.shift(0.3 * DOWN))
        subcode = Tex("is a subcode of (1)", font_size=26).next_to(
            code_list[2][1], DOWN, buff=0.15, aligned_edge=LEFT
        )
        self.play(Write(subcode))
        self.next_slide()

        self.play(FadeOut(code_list[2]), FadeOut(subcode))
        self.play(code_list[1].animate.shift(0.3 * DOWN))
        pushforward_note = VGroup(
            colored_tex(
                "Only valid if $", (r"\varphi", BLUE), "$ is injective on $",
                ("D", YELLOW), "$,",
                font_size=26,
            ),
            colored_tex("which forces $", ("n", YELLOW), "$ to be small.", font_size=26),
            colored_tex(
                "The resulting code has parameters $[", ("n", YELLOW), r",\leq\!",
                ("k", RED), "]$.",
                font_size=26,
            ),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        pushforward_note.next_to(code_list[1][1], DOWN, buff=0.15, aligned_edge=LEFT)
        self.play(Write(pushforward_note))
        self.next_slide()

        # One line, built from two Texes: as a single Tex, LaTeX wraps it at its
        # text width. The second is placed after the first and shifted so their
        # baselines match (the bottoms of "A" in Assuming and "S" in See).
        conditions = Tex(
            r"*Assuming $\ell\mid q-1$ and $\ker\varphi$"
            r" is cyclic and lies in $\mathcal{E}_1(\mathbb{F}_q)$.",
            font_size=18,
        )
        proofs = Tex("See the report for proofs.", font_size=18)
        proofs.next_to(conditions, RIGHT, buff=0.1)
        proofs.shift((conditions[0][1].get_bottom()[1] - proofs[0][0].get_bottom()[1]) * UP)
        assumptions = VGroup(conditions, proofs).to_edge(DOWN, buff=0.3).set_x(0)

        # Focus on (1), the pullback code, and its properties.
        self.play(FadeOut(code_list[1]), FadeOut(pushforward_note), Unwrite(action))
        self.play(code_list[0].animate.to_edge(UP, buff=0.5), Write(assumptions))

        # ell = #ker(phi) is blue like the isogeny.
        ell = (r"\ell", BLUE)

        def color_fraction(fraction, numerator_colors, denominator_colors):
            """Colour a \\tfrac's glyphs left to right, above and below its bar."""
            bar = max(fraction, key=lambda glyph: glyph.width)
            for side, colors in ((1, numerator_colors), (-1, denominator_colors)):
                glyphs = sorted(
                    (g for g in fraction if g is not bar and side * (g.get_y() - bar.get_y()) > 0),
                    key=lambda glyph: glyph.get_x(),
                )
                for glyph, color in zip(glyphs, colors):
                    glyph.set_color(color)

        def bullet(*lines):
            body = VGroup(*lines).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
            mark = MathTex(r"\bullet", font_size=26).next_to(body[0], LEFT, buff=0.2)
            return VGroup(mark, body)

        size = 26
        rate_line = colored_tex(
            "preserves the rate $", r"\tfrac{k}{n}", "$ and relative designed distance",
            font_size=size,
        )
        distance_line = colored_tex(
            "$", r"\tfrac{n-k}{n}", "$ of the original code", font_size=size
        )
        color_fraction(rate_line[1], [RED], [YELLOW])
        color_fraction(distance_line[1], [YELLOW, WHITE, RED], [YELLOW])
        # A part can't end inside a subscript, so C_ell is one part; its
        # rightmost glyph is the ell.
        matrix_product_line = colored_tex(
            r"is a matrix-product code $(C_1,\hdots,", r"C_{\ell}", r")\mathbf{V}$, where",
            font_size=size,
        )
        max(matrix_product_line[1], key=lambda glyph: glyph.get_x()).set_color(BLUE)
        properties = VGroup(
            bullet(colored_tex(
                "has parameters $[", ell, ("n", YELLOW), ",", ell, ("k", RED),
                r",\geq\!", ell, "(", ("n", YELLOW), "-", ("k", RED), ")]$, where $",
                ell, r"=\#\ker", (r"\varphi", BLUE), "$",
                font_size=size,
            )),
            bullet(rate_line, distance_line),
            bullet(
                colored_tex("is quasi-cyclic of co-index $", ell, "$", font_size=size),
                colored_tex(
                    "(i.e.\\ invariant under $", ell, "$-fold cyclic shifts)",
                    font_size=size,
                ),
            ),
            bullet(
                matrix_product_line,
                Tex(r"each $C_j$ is a `twist' of the original elliptic code", font_size=size),
                Tex(r"and $\mathbf{V}$ is a Vandermonde matrix", font_size=size),
            ),
        ).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        properties.next_to(code_list[0], DOWN, buff=0.5).set_x(0)

        for item in properties:
            self.play(Write(item))
            self.next_slide()

        # Closing slide, after TeacherStudentsScene from 3b1b's manim: three
        # students (scaled down) in the bottom-left corner, the teacher flipped
        # to face them in the bottom-right corner.
        students = VGroup(*(
            pi_creature(mode, color, height=1.5)
            for mode, color in zip(("happy", "hooray", "connving"), (BLUE_D, BLUE_E, BLUE_C))
        )).arrange(RIGHT, buff=0.25).to_corner(DL, buff=0.4)
        teacher = pi_creature("wave_1", GREY_BROWN, height=1.9).flip().to_corner(DR, buff=0.4)
        thanks = Tex(r"\textbf{Thank you!}", font_size=72).move_to(1.4 * UP)
        self.wipe(self.mobjects_without_canvas, Group(thanks, students, teacher), direction=UP)
        self.wait(0.1)