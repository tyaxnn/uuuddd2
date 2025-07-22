import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
import os
import sys
from datetime import datetime

def height_map_plot(kx,ky,data,discription,plot_title,output_dir):
    # 元のデータの範囲に基づいて、100x100のグリッドを定義
    grid_kx = np.linspace(kx.min(), kx.max(), 200)
    grid_ky = np.linspace(ky.min(), ky.max(), 200)
    grid_kx, grid_ky = np.meshgrid(grid_kx, grid_ky)

    grid_berry = griddata((kx, ky), data, (grid_kx, grid_ky), method='cubic')

    fig, ax = plt.subplots(figsize=(8, 7))

    # 塗りつぶし等高線プロット (色)
    contourf = ax.contourf(grid_kx, grid_ky, grid_berry, levels=10, cmap='viridis')

    # 等高線 (線)
    ax.contour(grid_kx, grid_ky, grid_berry, levels=contourf.levels, colors='white', linewidths=0.5, alpha=0.7)
    
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
    print(f"save as '{output_filename_num_band}'")
    plt.close()

def plot_berry_contour(file_path,size):

    if size != 2 and size != 6:
        print(f"invalid size")

    #------------------------------------------------------------
    #                       出力の準備
    #------------------------------------------------------------

    base_name = os.path.basename(file_path)
    name_without_ext, _ = os.path.splitext(base_name)

    now = datetime.now()
    time_str = now.strftime("%Y%m%d%H%M")

    folder_name = time_str + "_" + name_without_ext

    os.makedirs("output/" + folder_name, exist_ok=True)
    output_dir = './output/' + folder_name + "/" 

    print(f"input : {file_path}")
    print(f"output: {output_dir}")

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                plot_title = lines[1].strip() # strip()で改行文字などを除去
    except Exception as e:
        print(f"警告: ファイル '{file_path}' からタイトル行を読み取れませんでした: {e}")

    try:
        data = np.loadtxt(file_path,skiprows=2, comments='#', delimiter=",")
        if data.shape[1] != 2 + size * 5:
            print(f"エラー: ファイル '{file_path}' は32列である必要があります。")
            return
    except Exception as e:
        print(f"エラー: ファイル '{file_path}' の読み込み中にエラーが発生しました: {e}")

    kx = data[:, 0]
    ky = data[:, 1]

    #------------------------------------------------------------
    #                       ベリー曲率のplot
    #------------------------------------------------------------

    for i in range(2,2 + size):

        berry = data[:, i]

        discription = "b" + str(i - 1)
        
        height_map_plot(kx,ky,berry,discription,plot_title,output_dir)

    #------------------------------------------------------------
    #                       ベリー曲率のplot
    #------------------------------------------------------------

    for i in range(2 + size,2 + size * 2):

        energy = data[:, i]

        discription = "e" + str(i - 1 - size)
        
        height_map_plot(kx,ky,energy,discription,plot_title,output_dir)

    #------------------------------------------------------------
    #                       spin期待値のplot
    #------------------------------------------------------------

    for i in range(2 + size * 2,2 + size * 3):

        spin_num = data[:, i]

        fig, ax = plt.subplots(figsize=(8, 7))

        ax.scatter(kx[spin_num==0], ky[spin_num==0], color='blue', label='Spin 0', s=0.5)

        ax.scatter(kx[spin_num==1], ky[spin_num==1], color='red', label='Spin 1', s=0.5)

        discription = "s" + str(i - 1 - size * 2)

        ax.set_title(f"{plot_title} {discription}")
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='--', alpha=0.5)

        output_filename_num_band = output_dir + discription + '.png'

        plt.savefig(output_filename_num_band, dpi=300, bbox_inches='tight')
        print(f"save as '{output_filename_num_band}'")
        plt.close()

    #------------------------------------------------------------
    #                       BCD xのplot
    #------------------------------------------------------------

    for i in range(2 + size * 3,2 + size * 5,2):

        bcdx = data[:, i]

        discription = "dx" + str(int((i - size * 3)/2))

        height_map_plot(kx,ky,bcdx,discription,plot_title,output_dir)

    #------------------------------------------------------------
    #                       BCD yのplot
    #------------------------------------------------------------

    for i in range(2 + size * 3,2 + size * 5,2):

        bcdx = data[:, i + 1]

        discription = "dy" + str(int((i - size * 3)/2))

        height_map_plot(kx,ky,bcdx,discription,plot_title,output_dir)


if __name__ == '__main__':
    # 'berry.dat'を読み込んでプロットを作成
    plot_berry_contour('dats/a2_500_tmduuuddd_lambda0p1_j0p25_mu-0p1.csv',6)