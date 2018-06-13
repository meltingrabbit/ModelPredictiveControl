function dx = RungeKuttaSS(x, u, t, dt, A, B)
	% 2018-06-12
	% ルンゲクッタ（状態空間モデル用）
	% 1ステップ回す
	% dist			: 外乱ベクトル


	k1 = A*x                + B*u;
	k2 = A*(x + k1*0.5*dt ) + B*u;
	k3 = A*(x + k2*0.5*dt ) + B*u;
	k4 = A*(x + k3*dt )     + B*u;

	dx = dt .* (k1 + 2*k2 + 2*k3 + k4) ./ 6.0;
end


% function dx = RungeKuttaSS(x, u, t, dt, A, B, dist)
% 	% 2018-06-12
% 	% ルンゲクッタ（状態空間モデル用）
% 	% 1ステップ回す
% 	% dist			: 外乱ベクトル


% 	k1 = A*x                + B*u + dist;
% 	k2 = A*(x + k1*0.5*dt ) + B*u + dist;
% 	k3 = A*(x + k2*0.5*dt ) + B*u + dist;
% 	k4 = A*(x + k3*dt )     + B*u + dist;

% 	dx = dt .* (k1 + 2*k2 + 2*k3 + k4) ./ 6.0;
% end



%{
void RungeKutta(vector<double (*)(double t, vector<double>& x)>& f, vector<double>& x, double t0, double tn, double dt) {
	int num = x.size();
	vector<double> k1(num), k2(num), k3(num), k4(num);		// RungeKuttaの係数
	vector<double> temp(num);
	double t = t0;

	// int count = 0;

	while (t < tn) {
		for (int j=0; j<num; ++j) {
			k1[j] = (*f[j])(t, x);
			temp[j] = x[j] + dt*k1[j]/2;
		}
		for (int j=0; j<num; ++j) {
			k2[j] = (*f[j])(t+dt/2, temp);
		}
		for (int j=0; j<num; ++j) {
			temp[j] = x[j] + dt*k2[j]/2;
		}
		for (int j=0; j<num; ++j) {
			k3[j] = (*f[j])(t+dt/2, temp);
		}
		for (int j=0; j<num; ++j) {
			temp[j] = x[j] + dt*k3[j];
		}
		for (int j=0; j<num; ++j) {
			k4[j] = (*f[j])(t+dt, temp);
			x[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*dt/6;
		}
		t += dt;
		// count++;
	}
	// pn(count);
}
%}

