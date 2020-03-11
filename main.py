#!/usr/bin/env python

# %%
import numpy as np
from scipy.integrate import solve_ivp
from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider, Range1d
from bokeh.plotting import (
    ColumnDataSource,
    figure,
    output_file,
    show,
    save,
    curdoc,
)
from bokeh.io.doc import set_curdoc

# %%
def seir_ode(t, y, beta, gamma, sigma, mu, nu):
    # y[0] : susceptible (S)
    # y[1] : exposed (E)
    # y[2] : infected (I)
    # y[3] : resistant (R)

    ntot = y[0] + y[1] + y[2] + y[3]
    array_return = [
        mu * (ntot - y[0]) - beta * y[0] * y[2] / ntot - nu * y[0],
        beta * y[0] * y[2] / ntot - (mu + sigma) * y[1],
        sigma * y[1] - (mu + gamma) * y[2],
        gamma * y[2] - mu * y[3] + nu * y[0],
    ]
    return array_return


# %%
def solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span):

    n_eval = 101

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)

    sol = solve_ivp(
        seir_ode, t_span, y0, t_eval=t_eval, args=(beta, gamma, sigma, mu, nu)
    )

    return sol


def plot_with_bokeh():

    beta = 0.8
    gamma = 0.1
    sigma = 0.5
    mu = 0.0
    nu = 0.0
    y0 = [10, 1, 0, 0]
    t_span = [0, 30]  # days

    sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)

    x = sol.t
    y_s, y_e, y_i, y_r = sol.y

    source = ColumnDataSource(data=dict(x=x, y_s=y_s, y_e=y_e, y_i=y_i, y_r=y_r,))

    p = figure(
        plot_width=800, plot_height=600, x_axis_label="Days", y_axis_label="Population",
    )

    p.line(
        "x",
        "y_s",
        source=source,
        line_width=3,
        color="orange",
        legend_label="Susceptible",
    )
    p.line(
        "x",
        "y_e",
        source=source,
        line_width=3,
        color="dodgerblue",
        legend_label="Exposed",
    )
    p.line(
        "x",
        "y_i",
        source=source,
        line_width=3,
        color="orangered",
        legend_label="Infected",
    )
    p.line(
        "x",
        "y_r",
        source=source,
        line_width=3,
        color="seagreen",
        legend_label="Resistant",
    )

    slider_beta = Slider(start=0.0, end=1, value=0.8, step=0.1, title="\u03B2",)
    slider_gamma = Slider(start=0.0, end=1, value=0.1, step=0.1, title="\u03B3")
    slider_sigma = Slider(start=0.0, end=1, value=0.5, step=0.1, title="\u03C3")
    slider_mu = Slider(start=0.0, end=1, value=0.0, step=0.1, title="\u03BC")
    slider_nu = Slider(start=0.0, end=1, value=0.0, step=0.1, title="\u03BD")

    slider_s = Slider(start=0, end=100, value=10, step=1, title="N(Susceptible)")
    slider_e = Slider(start=0, end=100, value=1, step=1, title="N(Exposed)")
    slider_i = Slider(start=0, end=100, value=0, step=1, title="N(Infected)")
    slider_r = Slider(start=0, end=100, value=0, step=1, title="N(Recovered)")
    slider_t = Slider(start=0, end=100, value=30, step=1, title="Duration (days)")

    def callback_beta(attr, old, new):
        sol = solve_seir_model(new, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_gamma(attr, old, new):
        sol = solve_seir_model(beta, new, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_sigma(attr, old, new):
        sol = solve_seir_model(beta, gamma, new, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_mu(attr, old, new):
        sol = solve_seir_model(beta, gamma, sigma, new, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_nu(attr, old, new):
        sol = solve_seir_model(beta, gamma, sigma, mu, new, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_s(attr, old, new):
        y0[0] = new
        sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_e(attr, old, new):
        y0[1] = new
        sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_i(attr, old, new):
        y0[2] = new
        sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_r(attr, old, new):
        y0[3] = new
        sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    def callback_t(attr, old, new):
        t_span[1] = new
        sol = solve_seir_model(beta, gamma, sigma, mu, nu, y0, t_span)
        source.data["x"] = sol.t
        source.data["y_s"] = sol.y[0]
        source.data["y_e"] = sol.y[1]
        source.data["y_i"] = sol.y[2]
        source.data["y_r"] = sol.y[3]

    slider_beta.on_change("value", callback_beta)
    slider_gamma.on_change("value", callback_gamma)
    slider_sigma.on_change("value", callback_sigma)
    slider_mu.on_change("value", callback_mu)
    slider_nu.on_change("value", callback_nu)

    slider_s.on_change("value", callback_s)
    slider_e.on_change("value", callback_e)
    slider_i.on_change("value", callback_i)
    slider_r.on_change("value", callback_r)

    slider_t.on_change("value", callback_t)

    # draw_plot()

    sliders_params = column(
        slider_beta, slider_gamma, slider_sigma, slider_mu, slider_nu
    )
    sliders_inits = column(slider_s, slider_e, slider_i, slider_r, slider_t)
    layout = column(p, row(sliders_params, sliders_inits),)

    curdoc().add_root(layout)


plot_with_bokeh()
