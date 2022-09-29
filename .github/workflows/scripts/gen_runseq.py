#!/usr/bin/env python3

class RunSeq:
    def __init__(self, outfile, mode='w'):
        self.__time_loop = list()
        self.__outfile_name = outfile
        self.__outfile = None
        self.__mode = mode

    def __enter__(self):
        self.__outfile = open(self.__outfile_name, self.__mode, encoding="utf-8")
        self.__outfile.write("runSeq:: \n")
        return self

    def __exit__(self, *_):
        self.__exit_sequence()
        self.__outfile.close()
        return False

    @property
    def time_loop(self):
        if self.__time_loop:
            return self.__time_loop[-1][0]  # this is the first element of the last tuple
        else:
            return 0

    @property
    def active_depth(self):
        if self.__time_loop:
            return self.__time_loop[-1][1]
        else:
            return -1

    def enter_time_loop(self, coupling_time, active=True, newtime=True, addextra_atsign=False):
        if newtime:
            if addextra_atsign:
                self.__outfile.write ("@@" + str(coupling_time) + " \n" )
            else:
                self.__outfile.write ("@" + str(coupling_time) + " \n" )
            if active:
                self.__time_loop.append((self.time_loop+1, self.active_depth+1))
            else:
                self.__time_loop.append((self.time_loop+1, self.active_depth))

    def add_action(self, action, if_add):
        if if_add:
            self.__outfile.write ("  {}\n".format(action))

    def leave_time_loop(self, leave_time, if_write_hist_rest=False, addextra_atsign=False ):
        if leave_time and self.__time_loop:
            _, active_depth = self.__time_loop.pop()
            #if if_write_hist_rest or active_depth == 0:
            #    self.__outfile.write ("  MED med_phases_history_write        \n" )
            #    self.__outfile.write ("  MED med_phases_restart_write        \n" )
            #    self.__outfile.write ("  MED med_phases_profile              \n" )
            if addextra_atsign: 
                self.__outfile.write ("@@ \n" )
            else:
                self.__outfile.write ("@ \n" )

    def __exit_sequence(self):
        while self.__time_loop:
            self.leave_time_loop(True)
        self.__outfile.write ("::\n" )
