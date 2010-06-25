%%% File: Fdtd.erl
%%% Description: FDTD 1-Dimentional Code

-module(fdtd).
%-export([main/0]).
-compile(export_all).

%% Physical Constant
-record(physics, {pi, light, permittivity, permeability}).

%% Analysis Area Info
-record(area, {px, nx, dx}).

%% Time Info.
-record(time, {cfl, stop_time, dt, last_step}).

calc_const() ->
		PI = 3.14159,
		LIGHT = 2.99792458e8,
		PERMITTIVITY = 8.85418782e-12,
		PERMEABILITY = 4.0 * PI * 1.0e-7,
		Physics = #physics{pi = PI, light = LIGHT,
											 permittivity = PERMITTIVITY,
											 permeability = PERMEABILITY},
		PX = 1.0,
		NX = 1000,
		DX = PX / NX,
		Area = #area{px = PX, nx = NX, dx = DX},
		CFL = 0.5,
		STOP_TIME = 30.0e-9,
		DT = CFL / (LIGHT * (1.0 / DX)),
		LAST_STEP = erlang:trunc(STOP_TIME / DT),
		Time = #time{cfl = CFL, stop_time = STOP_TIME, dt = DT, last_step = LAST_STEP},
		{Physics, Area, Time}.

print_analysis_info({_, Area, Time}) ->
		io:format("Area  X: ~f [m]~n", [Area#area.px]),
		io:format("Grid  X: ~B~n", [Area#area.nx]),
		io:format("Delta X: ~f [m]~n", [Area#area.dx]),
		io:format("CFL    : ~f\n", [Time#time.cfl]),
		io:format("dt     : ~e [sec]~n", [Time#time.dt]),
		io:format("Step   : ~B~n", [Time#time.last_step]).

create_list(N) ->
		[0.0 || _X <- lists:seq(1, N)].

set_soft_source(0, [H | T], Value) ->
		[H + Value | T];
set_soft_source(Index, [H | T], Value) ->
		[H | set_soft_source(Index - 1, T, Value)].

calc_stimulus(Step, Dt) ->
		Tc = 10.0e-9,
		Pw = 0.1e-9,
		T_current = (erlang:float(Step) - 0.5) * Dt,
		T_eff = (T_current - Tc) / Pw,
		1.0 * math:exp(-1.0 * (T_eff * T_eff)).

create_diff_list(L) ->
		[_ | T] = L,
		[_ | T1] = lists:reverse(L),
		NewL = lists:reverse(T1),
		[X - Y || {X, Y} <- lists:zip(T, NewL)].		

cut_both_end_element([_|T]) ->
		[_ | T1] = lists:reverse(T),
		lists:reverse(T1).

calc_ey(Const, Ey, Hz) ->
		Diff_Hz = create_diff_list(Hz),
		Cut_Ey = cut_both_end_element(Ey),
		[send_calc_ey_message({E, DH}, Const) || {E, DH} <- lists:zip(Cut_Ey, Diff_Hz)].

send_calc_ey_message(EH, Const) ->
		Pid = spawn(fdtd, calc_ey_main, []),
		Pid ! {self(), {EH, Const}},
		receive
				{Pid, Msg} -> Msg
		end.

calc_ey_main() ->
		receive
				{From, {{E, DH}, {Eps, Dt, Dx}}} ->
						From ! {self(), E - Dt / (Eps * Dx) * DH}
		end.

calc_hz(Const, Hz, Ey) ->
		Diff_Ey = create_diff_list(Ey),
		[send_calc_hz_message({H, DE}, Const) || {H, DE} <- lists:zip(Hz, Diff_Ey)].

send_calc_hz_message(EH, Const) ->
		Pid = spawn(fdtd, calc_hz_main, []),
		Pid ! {self(), {EH, Const}},
		receive
				{Pid, Msg} -> Msg
		end.
						
calc_hz_main() ->
		receive
				{From, {{H, DE}, {Mu, Dt, Dx}}} ->
						From ! {self(), H - Dt / (Mu * Dx) * DE}
		end.

ey_boundary({Light, Dt, Dx}, {Ey_at_1, Ey_at_N1}, Hist_Ey) ->
		C = Light * Dt,
		[Hist_Ey_0, Hist_Ey_1, Hist_Ey_2, Hist_Ey_3] = Hist_Ey,
		Ey_0 = Hist_Ey_1 + (C - Dx) / (C + Dx) * (Ey_at_1 - Hist_Ey_0),
		Ey_N = Hist_Ey_2 + (C - Dx) / (C + Dx) * (Ey_at_N1 - Hist_Ey_3),
		{Ey_0, Ey_N}.

calc(Const, EH_Field, Step) ->
		{Physics, Area, Time} = Const,
		{Ey, Hz, Hist_Ey} = EH_Field,
		Stimulus = calc_stimulus(Step, Time#time.dt),
		Hz1 = set_soft_source(500, Hz, Stimulus * Time#time.dt),
		Ey1 = calc_ey({Physics#physics.permittivity, Time#time.dt, Area#area.dx}, Ey, Hz1),
		{Ey_at_0, Ey_at_N} = ey_boundary({Physics#physics.light, Time#time.dt, Area#area.dx},
																		 {lists:nth(1, Ey1), lists:last(Ey1)}, Hist_Ey),
		Ey2 = [Ey_at_0] ++ Ey1 ++ [Ey_at_N],
		Hist_Ey1 = [lists:nth(1, Ey2), lists:nth(2, Ey2), 
								lists:nth(length(Ey2) - 1, Ey2), lists:nth(length(Ey2), Ey2)],
		Hz2 = calc_hz({Physics#physics.permeability, Time#time.dt, Area#area.dx}, Hz1, Ey2),
		{Ey2, Hz2, Hist_Ey1, Stimulus}.

print_point(Fd, Time, Stim, V1, V2) ->
		io:fwrite(Fd, "~g, ~g, ~g, ~g~n", [Time, Stim, V1, V2]).

print_line(Fd, Step, Hz) ->
		io:fwrite(Fd, "# Step = ~B~n", [Step]),
		[io:fwrite(Fd, "~B, ~g~n", [I, Data]) || 
				{I, Data} <- lists:zip(lists:seq(0, length(Hz)-1), Hz)],
		io:fwrite(Fd, "~n",[]).

fdtd_loop(Const, EH_Field, Step, {Fd1, Fd2}) ->
		{_, _, Time} = Const,
		Last_Step = Time#time.last_step,
		_ = if Step rem 500 == 0 -> io:format("Step: ~B~n", [Step]); true -> 0 end,
		case Step of
				Last_Step -> file:sync(Fd1), file:close(Fd1),
										 file:sync(Fd2), file:close(Fd2);
				_ -> {Ey, Hz, Hist_Ey, Stimulus} = calc(Const, EH_Field, Step),
						 print_point(Fd1, Step * Time#time.dt, Stimulus, 
												 lists:nth(500, Hz), lists:nth(900, Hz)),
						 _ = if Step rem 100 == 0 -> print_line(Fd2, Step, Hz); true -> 0 end,
						 fdtd_loop(Const, {Ey, Hz, Hist_Ey}, Step+1, {Fd1, Fd2})
 	  end.

calc_fdtd(Const) ->
		{_, Area, _} = Const,
		Ey = create_list(Area#area.nx + 1),
		Hz = create_list(Area#area.nx),
		Hist_Ey = create_list(4),
		{ok, Fd1} = file:open("fdtd_erlang_point.csv", [write]),
		{ok, Fd2} = file:open("fdtd_erlang_line.csv", [write]),
		fdtd_loop(Const, {Ey, Hz, Hist_Ey}, 0, {Fd1, Fd2}).

main() ->
		Const = calc_const(),
		print_analysis_info(Const),
		calc_fdtd(Const).


