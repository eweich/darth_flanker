from smile.common import *

import os
import random
import copy
from math import cos, sin, sqrt, radians, log, pi
from numpy import linspace

# FLANKER VARIABLES

num_trials = 1 # (len(evidence_conditions)-1) * 4 + 2 * num_trials
num_blocks = 1
mouse = False
conds = ['congruent', 'incongruent', 'neutral']
resp_keys = ['G', 'H']

# CHANGE THIS LINE FOR G, H
pos_resps = resp_keys

num_flanks = 2
evidence_conditions = [0., 45.]
num_locs = 8
def_sat=[255.,255.,255.]
config_df = 10.
line_width = 1.
num_reps = 1

FONT_SIZE = .05#6
INST_FONT_SIZE = .025
orient_font_size = .05
item_font_size = .015

# Around 2.3 seconds each
orient_dur = .4
wait_dur = .5
wait_jitter = .5
resp_delay = .2
iti = .75
iti_jitter = .5
feedback_dur = 1.0


# list generation
def gen_fblocks():
    flanker_blocks = []
    for blo in range(num_blocks):
        for b in range(num_reps):
            temp_block = []
            for i in range(num_trials):
                for l in range(num_locs):
                    for s in evidence_conditions:
                        if s == 0.:
                            # NEUTRAL 1 RIGHT
                            trial = {'condition': "=",
                                     'target_ev': 45.,
                                     'inner_flank_ev': 45.,
                                     'outer_flank_ev':-45.,
                                     'loc':(360./num_locs)*l,
                                     'corr_resp': pos_resps[-1]}
                            temp_block.append(trial.copy())


                            #Neutral 2 RIGHT
                            trial['inner_flank_ev'] = -45.
                            trial['outer_flank_ev'] = 45.
                            trial['condition'] = "~"
                            temp_block.append(trial.copy())

                            # Neutral 2 LEFt
                            trial['corr_resp'] = pos_resps[0]
                            trial['inner_flank_ev'] = 45.
                            trial['outer_flank_ev'] = -45.
                            trial['target_ev'] = -45.
                            temp_block.append(trial.copy())

                            # Neutral 1 LEFT
                            trial['inner_flank_ev'] = -45.
                            trial['outer_flank_ev'] = 45.
                            trial['condition'] = "="
                            temp_block.append(trial.copy())
                        else:
                            # CONGRUENT RIGHT
                            trial = {'condition': "+",
                                     'target_ev': 45.,
                                     'inner_flank_ev': s,
                                     'outer_flank_ev': s,
                                     'loc':(360./num_locs)*l,
                                     'corr_resp': pos_resps[-1]}
                            temp_block.append(trial.copy())

                            #INCONGRUENT RIGHT
                            trial['inner_flank_ev'] = -1*s
                            trial['outer_flank_ev'] = -1*s
                            trial['condition'] = "-"
                            temp_block.append(trial.copy())

                            # INCONGRUENT LEFt
                            trial['corr_resp'] = pos_resps[0]
                            trial['inner_flank_ev'] = s
                            trial['outer_flank_ev'] = s
                            trial['target_ev'] = -45.
                            temp_block.append(trial.copy())

                            # CONGRUENT LEFT
                            trial['inner_flank_ev'] = -1*s
                            trial['outer_flank_ev'] = -1*s
                            trial['condition'] = "+"
                            temp_block.append(trial.copy())
            random.shuffle(temp_block)
            flanker_blocks.append(temp_block)
    return flanker_blocks


@Subroutine
def Flanks_with_resp(self, target_ev, num_flanks, inner_flank_ev, outer_flank_ev,
                     df, line_width, center_x, center_y, screen_h,
                     screen_w, corr_resp, loc, sat):
    self.X0 = center_x + ((screen_w/5.)*Ref(cos, loc*(pi/180.)))
    self.Y0 = center_y + ((screen_h/5.)*Ref(sin, loc*(pi/180.)))
    self.Yp = center_y + config_df + ((screen_h/5.)*Ref(sin, loc*(pi/180.)))
    self.Yn = center_y - config_df + ((screen_h/5.)*Ref(sin, loc*(pi/180.)))
    self.Xp = center_x + config_df + ((screen_w/5.)*Ref(cos, loc*(pi/180.)))
    self.Xn = center_x - config_df + ((screen_w/5.)*Ref(cos, loc*(pi/180.)))
    #self.bottom = center_y - config_df - (sqrt(8)*config_df*Ref(sin, Ref(radians, 45.)))
    #self.top = center_y + config_df + (sqrt(8)*config_df*Ref(sin, Ref(radians, 45.)))

    #self.center_x = center_x
    #self.center_y = center_y
    self.target_alpha = Ref(radians, target_ev)
    self.inner_flank_alpha = Ref(radians, inner_flank_ev)
    self.outer_flank_alpha = Ref(radians, outer_flank_ev)
    self.ev = [self.outer_flank_alpha]+[self.inner_flank_alpha] + [self.target_alpha] + [self.inner_flank_alpha] + [self.outer_flank_alpha]
    self.vert_flanks = [self.outer_flank_alpha]+[self.inner_flank_alpha]
    self.inner_flanks = [self.outer_flank_alpha]
    #self.X0 = center_x
    # Numflanks is doubled.
    with Parallel() as flanker:
        with Loop(self.ev) as lp:
            self.Xi = (self.X0) - (6.*config_df*(lp.i - num_flanks))
            with If(lp.i==num_flanks):
                self.col=[255.,255.,255.]
            with Else():
                self.col=sat

            with flanker.insert():
                line1 = Line(points=[self.Xi - (sqrt(2)*config_df*Ref(cos,lp.current)), self.Yp + (sqrt(2)*config_df*Ref(sin,lp.current)),
                                     self.Xi + (sqrt(2)*config_df*Ref(cos,lp.current)), self.Yp - (sqrt(2)*config_df*Ref(sin,lp.current))],
                             cap="square", width=line_width,color=self.col)
                Line(points=[self.Xi - (sqrt(2)*config_df*Ref(cos,lp.current)), self.Yn - (sqrt(2)*config_df*Ref(sin,lp.current)),
                              self.Xi + (sqrt(2)*config_df*Ref(cos,lp.current)), self.Yn + (sqrt(2)*config_df*Ref(sin,lp.current))],
                      cap="square", width=line_width,color=self.col)

        with Loop(self.vert_flanks) as vf:
            self.Yi1 = (self.Y0) - (6.*config_df*(vf.i-num_flanks))
            self.Yi2 = (self.Y0) + (6.*config_df*(vf.i-num_flanks))
            with flanker.insert():
                Line(points=[self.Xp - (sqrt(2)*config_df*Ref(cos, vf.current))-config_df, self.Yi1 + (sqrt(2)*config_df*Ref(sin,vf.current))+config_df,
                             self.Xp + (sqrt(2)*config_df*Ref(cos, vf.current))-config_df, self.Yi1 - (sqrt(2)*config_df*Ref(sin,vf.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xn - (sqrt(2)*config_df*Ref(cos,vf.current))+config_df, self.Yi1 - (sqrt(2)*config_df*Ref(sin,vf.current))-config_df,
                             self.Xn + (sqrt(2)*config_df*Ref(cos,vf.current))+config_df, self.Yi1 + (sqrt(2)*config_df*Ref(sin,vf.current))-config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xp - (sqrt(2)*config_df*Ref(cos, vf.current))-config_df, self.Yi2 + (sqrt(2)*config_df*Ref(sin,vf.current))+config_df,
                             self.Xp + (sqrt(2)*config_df*Ref(cos, vf.current))-config_df, self.Yi2 - (sqrt(2)*config_df*Ref(sin,vf.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xn - (sqrt(2)*config_df*Ref(cos,vf.current))+config_df, self.Yi2 - (sqrt(2)*config_df*Ref(sin,vf.current))-config_df,
                             self.Xn + (sqrt(2)*config_df*Ref(cos,vf.current))+config_df, self.Yi2 + (sqrt(2)*config_df*Ref(sin,vf.current))-config_df],
                             cap="square", width=line_width, color = self.col)
        with Loop(self.inner_flanks) as ifl:
            self.Yin1 = (self.Y0) - (6.*config_df*(ifl.i - (num_flanks-1)))
            self.Yin2 = (self.Y0) + (6.*config_df*(ifl.i - (num_flanks-1)))
            self.Xin1 = (self.X0) - (6.*config_df*(ifl.i - (num_flanks-1)))
            self.Xin2 = (self.X0) + (6.*config_df*(ifl.i - (num_flanks-1)))

            with flanker.insert():
                Line(points=[self.Xin1 - (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin1 + (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df,
                             self.Xin1 + (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin1 - (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin1 - (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin1 - (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df,
                             self.Xin1 + (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin1 + (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin2 - (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin2 + (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df,
                             self.Xin2 + (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin2 - (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin2 - (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin2 - (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df,
                             self.Xin2 + (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin2 + (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin1 - (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin2 + (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df,
                             self.Xin1 + (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin2 - (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin1 - (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin2 - (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df,
                             self.Xin1 + (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin2 + (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin2 - (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin1 + (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df,
                             self.Xin2 + (sqrt(2)*config_df*Ref(cos, ifl.current)), self.Yin1 - (sqrt(2)*config_df*Ref(sin,ifl.current))+config_df],
                             cap="square", width=line_width, color = self.col)
                Line(points=[self.Xin2 - (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin1 - (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df,
                             self.Xin2 + (sqrt(2)*config_df*Ref(cos,ifl.current)), self.Yin1 + (sqrt(2)*config_df*Ref(sin,ifl.current))-config_df],
                             cap="square", width=line_width, color = self.col)
    # wait for a response
    with UntilDone():
        Wait(resp_delay)
        Wait(until=line1.appear_time)
        # set what the correct response is for the trial
        kp = KeyPress(keys=resp_keys,
                      base_time=line1.appear_time['time'],
                      correct_resp=corr_resp)
        self.resp = kp.pressed
        self.press_time = kp.press_time
        self.rt = kp.rt
        self.correct = kp.correct


def _get_score(corr_trials, num_trials, rt_trials):
    """iRT = [np.log(r+1) for r in con_rt]
    cRT = [np.log(r+1) for r in incon_rt]
    iacc = con_corr/num_con
    cacc = incon_corr/num_incon"""

    g = (((corr_trials/num_trials)-0.5)/0.5)*\
        ( ((1/(sum([log(r+1) for r in rt_trials])/len(rt_trials)))-(1/log(2.5))) / \
        ((1/log(1.4))-(1/log(2.5))) )

    return int(g*100)

@Subroutine
def Flanker(self, screen_h, screen_w):

    self.flanker_instruct_top_text = 'You will be presented with a group of symbols. \n' + \
                                ' You will be asked to indicate the direction that the \n' + \
                                ' [b]middle arrow is pointing[/b], while ignoring any other symbols. \n\n' + \
                                ' [b]Press any key to proceed.[/b]'
    self.flanker_instruct_1 = '[b]Practice 1:[/b] \n \n' + \
                              'Look at the arrows as they appear inside the red circle. \n' + \
                              'Press the '+str(pos_resps[0])+' key if the arrow is pointing [b]left[/b]. \n' + \
                              'Press the '+str(pos_resps[1])+' key if the arrow is pointing [b]right[/b].'
    self.flanker_instruct_2 = '[b]Practice 2:[/b] \n \n' + \
                              'Respond to the arrow in the red circle while ignoring the other symbols. \n' + \
                              'Press the '+str(pos_resps[0])+' key if the arrow is pointing [b]left[/b]. \n' + \
                              'Press the '+str(pos_resps[1])+' key if the arrow is pointing [b]right[/b].'
    self.flanker_instruct_3 = '[b]Practice 3:[/b] \n \n' + \
                              'Try responding to the center arrow without the red circle. \n' + \
                              'Press the '+str(pos_resps[0])+' key if the arrow is pointing [b]left[/b]. \n' + \
                              'Press the '+str(pos_resps[1])+' key if the arrow is pointing [b]right[/b].'
    self.flanker_instruct_4 = '[b]Practice 4:[/b] \n \n' + \
                              'Now, the group of symbols will appear in different locations around the screen. \n' + \
                              'Look at the arrow inside the red circle and indicate its direction. \n' + \
                              'Press the '+str(pos_resps[0])+' key if the arrow is pointing [b]left[/b]. \n' + \
                              'Press the '+str(pos_resps[1])+' key if the arrow is pointing [b]right[/b]. \n\n' + \
                              'Press any key to begin.'
    self.flanker_instruct_5 = '[b]Practice 5:[/b] \n \n' + \
                              'The final practice set will look like the actual task. \n' + \
                              'The group of symbols will appear at different locations around the screen. \n' + \
                              'Find the arrow at the center of the group, and indicate its direction. \n' + \
                              'In between responses, please direct your gaze to the cross at the center of the screen. \n' + \
                              'Press the '+str(pos_resps[0])+' key if the arrow is pointing [b]left[/b]. \n' + \
                              'Press the '+str(pos_resps[1])+' key if the arrow is pointing [b]right[/b]. \n\n' + \
                              'Press any key to begin.'

    with Parallel():
        Label(text=self.flanker_instruct_top_text,
          markup=True,
          halign='center',
          #bottom=screen_h*.7,
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)

        # upper left
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=45., inner_flank_ev=45.,
                                outer_flank_ev=45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) - screen_w/4.,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) + screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0., sat=[255.,255.,255.])
        # upper middle
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=45., inner_flank_ev=45.,
                                outer_flank_ev=-45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)),
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) + screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0., sat=[255.,255.,255.])
        # upper right
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=45., inner_flank_ev=-45.,
                                outer_flank_ev=-45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) + screen_w/4.,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) + screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0., sat=[255.,255.,255.])

        # lower left
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=-45., inner_flank_ev=-45.,
                                outer_flank_ev=-45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) - screen_w/4.,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) - screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0, sat=[255.,255.,255.])
        # lower middle
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=-45., inner_flank_ev=45.,
                                outer_flank_ev=-45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)),
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) - screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0, sat=[255.,255.,255.])
        # lower right
        Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=-45., inner_flank_ev=45.,
                                outer_flank_ev=45.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) + screen_w/4.,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) - screen_h/4.,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=None,
                                loc=0., sat=[255.,255.,255.])
    with UntilDone():
        KeyPress()


    with Loop([[45.,pos_resps[1]],[-45.,pos_resps[0]],
               [45.,pos_resps[1]],[45.,pos_resps[1]]]) as prac_ev:
        Wait(1.0)
        p1 = Flanks_with_resp(num_flanks=2,
                                target_ev=prac_ev.current[0], inner_flank_ev=0.,
                                outer_flank_ev=0.,
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) ,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) ,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=prac_ev.current[1],
                                loc=0., sat=[0.,0.,0.])
        with If(p1.correct):
            Label(text='Correct!',color='green',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
        with Else():
            Label(text='Incorrect',color='red',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
    with Meanwhile():
        with Parallel():
            Ellipse(color='red', size=(80,80))
            Ellipse(color='black', size=(70,70))
            Label(text=self.flanker_instruct_1,
                  markup=True,
                  halign='center',
                  center_x=screen_w/2.,
                  center_y=screen_h/2.+screen_h/4.,
                  font_size=screen_h*INST_FONT_SIZE)

    Label(text="Press any key to proceed to the next practice.",
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    Wait(1.0)

    with Loop([[45.,pos_resps[1],45.,45.],[-45.,pos_resps[0],-45.,45.],
               [45.,pos_resps[1],-45.,-45.],[45.,pos_resps[1],-45.,45.]]) as prac_ev:
        Wait(1.0)
        p2 = Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=prac_ev.current[0],
                                inner_flank_ev=prac_ev.current[2],
                                outer_flank_ev=prac_ev.current[3],
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) ,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) ,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=prac_ev.current[1],
                                loc=0., sat=[255.,255.,255.])
        with If(p2.correct):
            Label(text='Correct!',color='green',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
        with Else():
            Label(text='Incorrect',color='red',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
    with Meanwhile():
        with Parallel():
            Ellipse(color='red', size=(80,80))
            Ellipse(color='black', size=(70,70))
            Label(text=self.flanker_instruct_2,
                  markup=True,
                  halign='center',
                  center_x=screen_w/2.,
                  center_y=screen_h/2.+screen_h/4.,
                  font_size=screen_h*INST_FONT_SIZE)

    Label(text="Press any key to proceed to the next practice.",
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    Wait(1.0)
    with Loop([[-45.,pos_resps[0],45.,45.],[45.,pos_resps[1],-45.,45.],
               [-45.,pos_resps[0],-45.,-45.],[-45.,pos_resps[0],-45.,45.],
               [45.,pos_resps[1],45.,45.],[45.,pos_resps[1],-45.,-45.]]) as prac_ev:
        Wait(1.0)
        p3 = Flanks_with_resp(num_flanks=num_flanks,
                                target_ev=prac_ev.current[0],
                                inner_flank_ev=prac_ev.current[2],
                                outer_flank_ev=prac_ev.current[3],
                                df=config_df,
                                line_width=line_width,
                                center_x=screen_w/2. - ((screen_w/5.)*Ref(cos, 0.)) ,
                                center_y=screen_h/2. - ((screen_h/5.)*Ref(sin, 0.)) ,
                                screen_h=screen_h,
                                screen_w=screen_w,
                                corr_resp=prac_ev.current[1],
                                loc=0., sat=[255.,255.,255.])
        with If(p3.correct):
            Label(text='Correct!',color='green',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
        with Else():
            Label(text='Incorrect',color='red',duration=feedback_dur, center_y=screen_h/2.-(screen_h/8.),font_size=screen_h*FONT_SIZE)
    with Meanwhile():
        Label(text=self.flanker_instruct_3,
              markup=True,
              halign='center',
              center_x=screen_w/2.,
              center_y=screen_h/2.+screen_h/4.,
              font_size=screen_h*INST_FONT_SIZE)
    Label(text="Press any key to proceed to the next practice.",
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    Wait(1.0)
    Label(text=self.flanker_instruct_4,
          markup=True,
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.+screen_h/4.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    with Loop([[45.,pos_resps[1],45.,(360./num_locs)*5,-45.],[-45.,pos_resps[0],-45,(360./num_locs)*1,-45.],
               [-45.,pos_resps[0],45.,(360./num_locs)*3,45.],[45.,pos_resps[1],45.,(360./num_locs)*6,45.],
               [45.,pos_resps[1],-45.,(360./num_locs)*2,45.],[-45.,pos_resps[0],-45.,(360./num_locs)*4,-45.]]) as prac_ev:
        Wait(1.0)
        p4 = Flanks_with_resp(num_flanks=num_flanks,
                            target_ev=prac_ev.current[0],
                            inner_flank_ev=prac_ev.current[2],
                            outer_flank_ev=prac_ev.current[4],
                            df=config_df,
                            line_width=line_width,
                            center_x=screen_w/2.,
                            center_y=screen_h/2.,
                            screen_h=screen_h,
                            screen_w=screen_w,
                            corr_resp=prac_ev.current[1],
                            loc=prac_ev.current[3], sat=[255.,255.,255.])
        with Meanwhile():
            with Parallel():
                Ellipse(color='red', size=(80,80),
                        center_x=screen_w/2. + ((screen_w/5.)*Ref(cos, prac_ev.current[3]*(pi/180.))) ,
                        center_y=screen_h/2. + ((screen_h/5.)*Ref(sin, prac_ev.current[3]*(pi/180.))))
                Ellipse(color='black', size=(70,70),
                        center_x=screen_w/2. + ((screen_w/5.)*Ref(cos, prac_ev.current[3]*(pi/180.))) ,
                        center_y=screen_h/2. + ((screen_h/5.)*Ref(sin, prac_ev.current[3]*(pi/180.))))
        with If(p4.correct):
                Label(text='Correct!',color='green',duration=feedback_dur, center_y=screen_h/2.,font_size=screen_h*FONT_SIZE)
        with Else():
                Label(text='Incorrect',color='red',duration=feedback_dur, center_y=screen_h/2.,font_size=screen_h*FONT_SIZE)
    Label(text="Press any key to proceed to the next practice.",
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    Wait(1.0)
    Label(text=self.flanker_instruct_5,
          markup=True,
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.+screen_h/4.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    with Loop([[45.,pos_resps[1],-45.,(360./6)*6,45.],[45.,pos_resps[1],45,(360./6)*3,45.],
               [-45.,pos_resps[0],-45.,(360./6)*1,-45.],[-45.,pos_resps[0],45.,(360./6)*5,-45.],
               [45.,pos_resps[1],-45.,(360./6)*4,-45.],[-45.,pos_resps[0],45.,(360./6)*2,45.],
               [45.,pos_resps[1],-45.,(360./6)*3,-45.],[-45.,pos_resps[0],-45.,(360./6)*1,45.],
               [-45.,pos_resps[0],45.,(360./6)*4,-45.],[45.,pos_resps[1],45.,(360./6)*6,45.]]) as prac_ev:
        Wait(1.0)
        p5 = Flanks_with_resp(num_flanks=num_flanks,
                            target_ev=prac_ev.current[0],
                            inner_flank_ev=prac_ev.current[2],
                            outer_flank_ev=prac_ev.current[4],
                            df=config_df,
                            line_width=line_width,
                            center_x=screen_w/2.,
                            center_y=screen_h/2.,
                            screen_h=screen_h,
                            screen_w=screen_w,
                            corr_resp=prac_ev.current[1],
                            loc=prac_ev.current[3], sat=[255.,255.,255.])

        with If(p5.correct):
            Label(text='Correct!',color='green',duration=feedback_dur, center_y=screen_h/2.+screen_h/8.,font_size=screen_h*FONT_SIZE)

        with Else():
            Label(text='Incorrect',color='red',duration=feedback_dur, center_y=screen_h/2.+screen_h/8.,font_size=screen_h*FONT_SIZE)

    with Meanwhile():
        Label(text="+", font_size=screen_h*orient_font_size)
    Label(text="You've completed the practice!\nPress any key to proceed to the task.",
          halign='center',
          center_x=screen_w/2.,
          center_y=screen_h/2.,
          font_size=screen_h*INST_FONT_SIZE)
    with UntilDone():
        KeyPress()
    res = Func(gen_fblocks)
    self.blocks = res.result

    with Loop(self.blocks) as block:
        Wait(2.0)

        self.trials_corr = 0.
        self.trials_num = 0.
        self.trials_rt = []

        # present fixation cross

        # loop over trials within block
        with Loop(block.current) as trial:

            # wait some jittered amount
            Wait(duration=iti, jitter=iti_jitter)

            # present stimulus
            test_stim = Flanks_with_resp(num_flanks=num_flanks,
                               target_ev=trial.current['target_ev'],
                               inner_flank_ev=trial.current['inner_flank_ev'],
                               outer_flank_ev=trial.current['outer_flank_ev'],
                               df=self.df,
                               line_width=line_width, center_x=screen_w/2.,
                               center_y=screen_h/2., screen_h=screen_h, screen_w=screen_w,
                               corr_resp=trial.current['corr_resp'],
                               loc=trial.current['loc'],
                               sat=def_sat)

            with If(test_stim.correct):
                self.trials_corr = self.trials_corr + 1.
            self.trials_num = self.trials_num + 1.
            self.trials_rt = self.trials_rt + [test_stim.rt]

            # Log it
            Log(trial.current,
                name="FL",
                block=block.i,
                resp=test_stim.resp,
                press_time=test_stim.press_time,
                rt=test_stim.rt,
                #orient_time=orient.appear_time,
                correct=test_stim.correct)
        with Meanwhile():
            orient = Label(text="+", #duration=orient_dur,
                           font_size=screen_h*orient_font_size)




        Wait(1.0)
        self.block_score = Ref(_get_score, self.trials_corr, self.trials_num,
                                           self.trials_rt)

        with Parallel():
            Rectangle(width="700sp", height="500sp",
                      color=[144./255., 175./255., 197./255.])
            pbfbC = Label(text=Ref(str, self.block_score)+" Points!",
                          font_size=screen_h*FONT_SIZE)
            Label(text="Your Score for this block:",
                  font_size=screen_h*FONT_SIZE, bottom=pbfbC.top + 10.)
            Label(text="Press any key to continue!",
                  font_size=screen_h*FONT_SIZE, top=pbfbC.bottom - 10.)
        with UntilDone():
            KeyPress()
        Log(name="flanker_block_score",
            block_score=self.block_score)
        # wait before suddenly presenting next block
        Wait(wait_dur, jitter=wait_jitter)




exp = Experiment()
Label(text='',duration=1.)
Wait(1.)
Flanker(screen_w=exp.screen.width,
        screen_h=exp.screen.height)

exp.run()
