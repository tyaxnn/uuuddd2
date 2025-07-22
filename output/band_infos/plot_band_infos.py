import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
import os
import sys
from datetime import datetime

def plot_berry_contour(file_path):

    # --- ステップ1: 入力パスから動的に出力ファイル名を生成 ---
    # 'dats/berry_fukui_...dat' -> 'berry_fukui_...dat'
    base_name = os.path.basename(file_path)
    # 'berry_fukui_...dat' -> 'berry_fukui_...'
    name_without_ext, _ = os.path.splitext(base_name)

    now = datetime.now()
    # strftimeメソッドで指定した書式に変換
    time_str = now.strftime("%Y%m%d%H%M")

    folder_name = time_str + "_" + name_without_ext

    os.makedirs("output/" + folder_name, exist_ok=True)
    # 'berry_fukui_...' -> 'berry_fukui_....png'
    output_dir = './output/' + folder_name + "/" 

    print(f"入力ファイル: {file_path}")
    print(f"出力ファイル: {output_dir}")

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            # 2行目からパラメータを抽出し、タイトルを生成
            if len(lines) > 1:
                plot_title = lines[1].strip() # strip()で改行文字などを除去
    except Exception as e:
        print(f"警告: ファイル '{file_path}' からタイトル行を読み取れませんでした: {e}")

    # --- ステップ1: データの読み込み ---
    try:
        # スペースまたはカンマ区切りデータを想定
        data = np.loadtxt(file_path,skiprows=2, comments='#', delimiter=",")
        if data.shape[1] != 32:
            print(f"エラー: ファイル '{file_path}' は32列である必要があります。")
            return
    except Exception as e:
        print(f"エラー: ファイル '{file_path}' の読み込み中にエラーが発生しました: {e}")

    kx = data[:, 0]
    ky = data[:, 1]

    for i in range(2,8):

        berry = data[:, i]

        # --- ステップ2: 補間用のグリッドを作成 ---
        # 元のデータの範囲に基づいて、100x100のグリッドを定義
        grid_kx = np.linspace(kx.min(), kx.max(), 200)
        grid_ky = np.linspace(ky.min(), ky.max(), 200)
        grid_kx, grid_ky = np.meshgrid(grid_kx, grid_ky)

        # --- ステップ3: 散乱データをグリッド上に補間 ---
        # 'cubic'は滑らかな補間、'linear'は線形補間
        grid_berry = griddata((kx, ky), berry, (grid_kx, grid_ky), method='cubic')

        # --- ステップ4: 等高線プロットの作成 ---
        fig, ax = plt.subplots(figsize=(8, 7))

        # 塗りつぶし等高線プロット (色)
        contourf = ax.contourf(grid_kx, grid_ky, grid_berry, levels=10, cmap='viridis')

        # contour_levels = np.linspace(0.000,0.00027,8)

        # 等高線 (線)
        ax.contour(grid_kx, grid_ky, grid_berry, levels=contourf.levels, colors='white', linewidths=0.5, alpha=0.7)

        discription = "b" + str(i - 1)
        
        # --- ステップ5: ラベルとカラーバーの設定 ---
        ax.set_title(plot_title + " " + discription)
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_aspect('equal', adjustable='box') # アスペクト比を1:1に
        fig.colorbar(contourf, ax=ax, label='Berry Value ($\Omega$)')
        ax.grid(True, linestyle='--', alpha=0.5)

        output_filename_num_band = output_dir + discription + '.png'
        
        # --- ステップ6: プロットをファイルに保存 ---
        plt.savefig(output_filename_num_band, dpi=300, bbox_inches='tight')
        print(f"等高線プロットを '{output_filename_num_band}' として保存しました。")
        plt.close()

    for i in range(8,14):

        energy = data[:, i]

        # --- ステップ2: 補間用のグリッドを作成 ---
        # 元のデータの範囲に基づいて、100x100のグリッドを定義
        grid_kx = np.linspace(kx.min(), kx.max(), 200)
        grid_ky = np.linspace(ky.min(), ky.max(), 200)
        grid_kx, grid_ky = np.meshgrid(grid_kx, grid_ky)

        # --- ステップ3: 散乱データをグリッド上に補間 ---
        # 'cubic'は滑らかな補間、'linear'は線形補間
        grid_energy = griddata((kx, ky), energy, (grid_kx, grid_ky), method='cubic')

        # --- ステップ4: 等高線プロットの作成 ---
        fig, ax = plt.subplots(figsize=(8, 7))

        # 塗りつぶし等高線プロット (色)
        contourf = ax.contourf(grid_kx, grid_ky, grid_energy, levels=10, cmap='viridis')

        # contour_levels = np.linspace(0.000,0.00027,8)

        # 等高線 (線)
        ax.contour(grid_kx, grid_ky, grid_energy, levels=contourf.levels, colors='white', linewidths=0.5, alpha=0.7)

        discription = "e" + str(i - 7)
        
        # --- ステップ5: ラベルとカラーバーの設定 ---
        ax.set_title(plot_title + " " + discription)
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_aspect('equal', adjustable='box') # アスペクト比を1:1に
        fig.colorbar(contourf, ax=ax, label='Berry Value ($\Omega$)')
        ax.grid(True, linestyle='--', alpha=0.5)

        output_filename_num_band = output_dir + discription + '.png'
        
        # --- ステップ6: プロットをファイルに保存 ---
        plt.savefig(output_filename_num_band, dpi=300, bbox_inches='tight')
        print(f"等高線プロットを '{output_filename_num_band}' として保存しました。")
        plt.close()

    for i in range(14,20):

        spin_num = data[:, i]

        fig, ax = plt.subplots(figsize=(8, 7))

        # spin=0のデータを青色でプロット
        ax.scatter(kx[spin_num==0], ky[spin_num==0], color='blue', label='Spin 0', s=0.5)

        # spin=1のデータを赤色でプロット
        ax.scatter(kx[spin_num==1], ky[spin_num==1], color='red', label='Spin 1', s=0.5)

        discription = "s" + str(i - 13)

        ax.set_title(f"{plot_title} {discription}")
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='--', alpha=0.5)



        output_filename_num_band = output_dir + discription + '.png'

        # --- ステップ5: プロットをファイルに保存 ---
        plt.savefig(output_filename_num_band, dpi=300, bbox_inches='tight')
        print(f"プロットを '{output_filename_num_band}' として保存しました。")
        plt.close()

    for i in range(20,32,2):

        bcdx = data[:, i]
        bcdy = data[:, i + 1]

        # --- ステップ3: プロットの作成 ---
        fig, ax = plt.subplots(figsize=(8, 7))

        magnitudes = np.sqrt(bcdx**2 + bcdy**2)

        step = 17
        quiver_plot = ax.quiver(
            kx[::step], ky[::step], bcdx[::step], bcdy[::step],
            magnitudes[::step], # ★矢印の色をベクトルの大きさで変える
            cmap='jet',         # 色のマップ（viridis, plasmaなどもおすすめ）
            angles='xy',        # 矢印の角度を(x,y)成分から計算
            scale_units='xy',   # スケールをx,y軸の単位に合わせる
            scale=0.5,        # ★スケールをデータの値に合わせて小さく調整
            width=0.001
        )

        discription = "d" + str(int((i - 18)/2))

        # --- ステップ4: ラベルとカラーバーの設定 ---
        ax.set_title(f"{plot_title} {discription}")
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='--', alpha=0.5)

        output_filename_num_band = output_dir + discription + '.png'

        # --- ステップ5: プロットをファイルに保存 ---
        plt.savefig(output_filename_num_band, dpi=300, bbox_inches='tight')
        print(f"プロットを '{output_filename_num_band}' として保存しました。")
        plt.close()

if __name__ == '__main__':
    # 'berry.dat'を読み込んでプロットを作成
    plot_berry_contour('dats/a1_500_tmduuudddum_lambda0p1_j0p25_mu-0p1.csv')