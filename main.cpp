//多結晶粒組織形成のプログラム
//粒番号は1-nm
//粒番号nmを液相

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の設定

#define ND 100					//差分計算における計算領域一辺の分割数
#define N 21						//考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)
#define INX 400					//描画window１辺(x方向)のピクセルサイズ
#define INY 400					//描画window１辺(y方向)のピクセルサイズ

	int nd=ND, ndm=ND-1;	//計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
	int nm=N-1, nmm=N-2;	//考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
	double PI=3.141592, time1;	//π、計算カウント数
	double RR=8.3145;						//ガス定数
	double ph[N][ND][ND], ph2[N][ND][ND];	//フェーズフィールド、フェーズフィールド補助配列
  double aij[N][N];//勾配エネルギー係数
  double wij[N][N];//ペナルティー項の係数
  double tij[N][N];//粒界の易動度
  double eij[N][N];//粒界移動の駆動力
  int m00h[N][ND][ND];//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
  int n00h[ND][ND];//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数

	void ini000();	//初期場の設定サブルーチン
	void graph_a();	//グラフ描画サブルーチン
	void datsave();	//データ保存サブルーチン
	void datin();		//データ入力サブルーチン

//******* メインプログラム ******************************************
int main(void)
{
	int i, j, k, l, ii, jj, kk, ll, it;	//整数
	int ip, im, jp, jm;									//整数
	int n1, n2, n3;											//整数
	int n00;		//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数
	//int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
	double time1max;				//計算カウント数の最大値（計算終了カウント）
	double delt, L, b1;			//時間�
	double M1;			//粒界の易動度
	double W1;			//ペナルティー項の係数
	double K1;			//勾配エネルギー係数
	double E1;			//粒界移動の駆動力
	double temp;		//温度
	double sum1, sum2, sum3;//各種の和の作業変数
	double pddtt;		//フェーズフィールドの時間変化率

  double gamma;	//粒界エネルギ密度
  double delta;	//粒界幅（差分ブロック数にて表現）
  double amobi;	//粒界の易動度
  double vm0;		//モル体積

//****** 計算条件および物質定数の設定 ****************************************

	printf("delt(5.0)=  "); scanf(" %lf",&delt);//時間刻みの入力
//	delt=5.0;

	temp=1000.0; 							//温度(K)
	L=2000.0;									//計算領域の一辺の長さ(nm)
	b1=L/(double)ND*1.0e-9; 	//差分プロック１辺の長さ(m)

	vm0=7.0e-6;								//モル体積
  gamma=0.5*vm0/RR/temp/b1;	//粒界エネルギ密度（0.5J/m^2）を無次元化
  delta=7.0;								//粒界幅（差分ブロック数にて表現）

	K1=8.0*delta*gamma/PI/PI;	//勾配エネルギー係数[式(4.40)]
	W1=4.0*gamma/delta;				//ペナルティー項の係数[式(4.40)]
  amobi=1.;	 M1=amobi*PI*PI/(8.0*delta);	//粒界の易動度[式(4.40)]
  E1=50.0/RR/temp;					//粒界移動の駆動力

	time1=0.;  time1max=10000001.;//計算カウント数の初期値と最大値

//*** 式(4.32) - 式(4.35)の配列（K,W,M,E）の設定 (MPF0.cppの場合の設定)*********
//	for(i=1;i<=nm;i++){
//		for(j=1;j<=nm;j++){
//			wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=0.0;  eij[i][j]=0.0;
//			if( (i==nm)||(j==nm) ){eij[i][j]=E1; tij[i][j]=M1;}
//			if(i>j){eij[i][j]=-eij[i][j];}
//			if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
//		}
//	}

//*** 式(4.36) - 式(4.39)の配列（K,W,M,E）の設定 **************************
	for(i=1;i<=nm;i++){
		for(j=1;j<=nm;j++){
			wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=M1;  eij[i][j]=0.0;
			if( (i==nm)||(j==nm) ){eij[i][j]=E1;}
			if(i>j){eij[i][j]=-eij[i][j];}
			if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
		}
	}

//*** 初期場の設定と描画Window表示 *****************************************
	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//フェーズフィールド、および補助配列の初期化
			}
		}
	}
	ini000();//初期場の設定
	//datin();//初期場の入力
 	gwinsize(INX,INY); ginit(2); gsetorg(0,0);//描画Window表示

//**** シミュレーションスタート ******************************
start: ;

	if((((int)(time1) % 200)==0)) {datsave();}//一定繰返しカウント毎に場を保存
	//if(time1==200.) {datsave();}//特定の時間の場を保存
	if((((int)(time1) % 20)==0)) {graph_a();} //一定繰返しカウント毎に場を表示

//**** 各差分プロックにおけるn00h[i][j]とm00h[n00][i][j]を調査 *********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

//--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
			n00=0;
			for(ii=1;ii<=nm;ii++){
				if( (ph[ii][i][j]>0.0)||
             ( (ph[ii][i][j]==0.0)&&(ph[ii][ip][j]>0.0)||
                                    (ph[ii][im][j]>0.0)||
                                    (ph[ii][i][jp]>0.0)||
                                    (ph[ii][i][jm]>0.0)   ) ){
            n00++; m00h[n00][i][j]=ii;
						//printf("%d  ", n00);
				}
			}
			n00h[i][j]=n00;
//--------------------------------------------------------------------------
		}
	}


//***** 発展方程式の計算 **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			for(n1=1; n1<=n00h[i][j]; n1++){
				ii=m00h[n1][i][j];  pddtt=0.0;
				for(n2=1; n2<=n00h[i][j]; n2++){
					jj=m00h[n2][i][j];  sum1=0.0;
					for(n3=1; n3<=n00h[i][j]; n3++){
           	kk=m00h[n3][i][j];
           	sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(ph[kk][ip][j]+ph[kk][im][j]
                                                +ph[kk][i][jp]+ph[kk][i][jm]-4.0*ph[kk][i][j])
              		+(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];//[式(4.31)の一部]
					}
          pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])
                  *(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[ii][i][j]*ph[jj][i][j]));
																									//フェーズフィールドの発展方程式[式(4.31)]
				}
				ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;		//フェーズフィールドの時間発展（陽解法）
				if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//フェーズフィールドの変域補正
				if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
			}
		}//j
	}//i

	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				ph[k][i][j]=ph2[k][i][j];
			}
		}
	}

//*** フェーズフィールドの規格化補正 ***********************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
		}
	}

//*********************************************************************
	if(keypress()){return 0;}//キー待ち状態
	time1=time1+1.;  if(time1<time1max){goto start;}//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}

//************ 初期場(フェーズフィールド)の設定サブルーチン *************
void ini000()
{
	int i, j, k, l, it;		//整数
	int ii, jj, kk;				//整数
	int ip, im, jp, jm;		//整数
	int x1, y1, x1h[10], y1h[10];//初期核の座標
	double sum1, t, r0, phi, r;
 	srand(3.0); // 乱数初期化
 	//srand(time(NULL)); // 乱数初期化

	//x1h[1]=0.2*nd;   y1h[1]=0.2*nd;		//初期核１の座標設定
	//x1h[2]=0.75*nd;  y1h[2]=0.4*nd;		//初期核２の座標設定
	//x1h[3]=0.5*nd;   y1h[3]=0.75*nd;	//初期核３の座標設定

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(ii=1;ii<=nm-1;ii++){ph[ii][i][j]=0.0;}
			ph[nm][i][j]=1.0;//nm番目のフェーズフィールドを１に初期化
		}
	}

	r0=5.0;
	for(ii=1;ii<=nm-1;ii++){
		//x1=x1h[ii]; y1=y1h[ii];
		x1=nd*DRND(1); y1=nd*DRND(1);//初期核の位置
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				r=sqrt( (double(i-x1))*(double(i-x1))+(double(j-y1))*(double(j-y1)) );
				if(r<=r0){ ph[ii][i][j]=1.0;  ph[nm][i][j]=0.0; } //初期核位置のフェーズフィールドを設定
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }//フェーズフィールドの規格化補正
		}
	}

}

//******* フェーズフィールドの描画サブルーチン ***************************************
void graph_a()
{
	int i, j, k, ii, jj;											//整数
	double col, col_R, col_G, col_B;					//色
	int ixmin=0, iymin=0, igx, igy, irad0;		//スクリーン座標系の設定
	int ixmax=INX, iymax=INY;									//描画Window範囲
	double x, xmax, xmin, y, ymax, ymin, rad0;//規格化座標系の設定

  //gcls();//画面クリア
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;		//描画領域（規格化されている）
	printf("time %f\n",time1);						//計算カウント数の表示
	rad0=1./nd/2.;              irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//座標計算				//スクリーン座標系に変換
			ii=i; jj=j;  if(i==nd){ii=0;} if(j==nd){jj=0;}//周期的境界条件
				col=0.; for(k=0;k<=nm;k++){ col+=ph[k][ii][jj]*ph[k][ii][jj]; }	//Σp^2を計算
				col_R=col_G=col_B=col;	//フェーズフィールドを明暗にて設定
				if(ph[nm][ii][jj]>0.5){col_R=0.; col_G=0.; col_B=1.;}	//nm番目の方位の粒を青にする
				gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));	//色設定
				grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);	//中塗り四角形の描画
		}
	}
	swapbuffers();//画面スワップ
}

//************ データ保存サブルーチン *******************************
void datsave()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k;		//整数

	stream = fopen("test.dat", "a");	//書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%e  \n", time1);	//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				fprintf(stream, "%e   ", ph[k][i][j]);//フェーズフィールドの保存
			}
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

//*********** データ入力サブルーチン **************************
void datin()
{
	FILE		*datin0;//ストリームのポインタ設定
	int 		i, j, k;//整数

	datin0 = fopen("test.ini", "r");//読み込み元のファイルをオープン
	fscanf(datin0, "%lf", &time1);	//計算カウント数の読み込み
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				fscanf(datin0, "%lf  ", &ph[k][i][j]);//フェーズフィールドの読み込み
			}
		}
	}
	fclose(datin0);//ファイルをクローズ

}
