import os
from pathlib import Path
import geopandas as gpd
from shapely.geometry import box
import rasterio
import numpy as np
from tqdm import tqdm

INPUT_DSM_TILES_DIRS = [
    Path(r"path"),
    Path(r"path")
]

# 格网文件路径
GRID_SHP_PATH = Path(r"path\grid_tiles.shp")
# 输出镶嵌后的大瓦片文件夹
OUTPUT_MOSAIC_DIR = Path(r"output")

TILE_INDEX_PATH = OUTPUT_MOSAIC_DIR.parent / "tile_spatial_index_combined.gpkg"

# 羽化带宽占瓦片宽度的百分比
FEATHER_WIDTH_PERCENT = 0.15
# 为每个大瓦片增加的重叠区（单位：像素）
OVERLAP_IN_PIXELS = 64

OUTPUT_MOSAIC_DIR.mkdir(parents=True, exist_ok=True)


def create_feather_weights(height, width, feather_percent=0.1):
    """为单个瓦片创建羽化权重图。"""
    if feather_percent <= 0:
        return np.ones((height, width), dtype=np.float32)
    feather_h = int(height * feather_percent)
    feather_w = int(width * feather_percent)
    feather_h = max(1, min(feather_h, height // 2))
    feather_w = max(1, min(feather_w, width // 2))
    weight_map = np.ones((height, width), dtype=np.float32)
    for i in range(feather_h):
        alpha = (i + 1) / (feather_h + 1.0)
        weight_map[i, :] *= alpha
        weight_map[height - 1 - i, :] *= alpha
    for j in range(feather_w):
        alpha = (j + 1) / (feather_w + 1.0)
        weight_map[:, j] *= alpha
        weight_map[:, width - 1 - j] *= alpha
    return weight_map


def mosaic_tiles_for_grid(tile_files, output_path, output_bounds, feather_percent=0.1):
    """将指定的瓦片列表，按照给定的输出边界，用羽化方式镶嵌起来。"""
    if not tile_files:
        print(f"没有提供用于镶嵌的瓦片文件。")
        return False
    print(f"\n开始为 {output_path.name} 镶嵌 {len(tile_files)} 个瓦片...")
    first_tile_meta = None
    tile_info_list = []

    for tile_path in tqdm(tile_files, desc="  读取元数据"):
        try:
            with rasterio.open(tile_path) as src:
                if first_tile_meta is None:
                    first_tile_meta = src.meta.copy()
                if not np.allclose([abs(first_tile_meta['transform'].a), abs(first_tile_meta['transform'].e)],
                                   [abs(src.transform.a), abs(src.transform.e)]):
                    print(f"    警告: 分辨率不匹配，跳过 {Path(tile_path).name}")
                    continue
                tile_info_list.append(
                    {"path": str(tile_path), "bounds": src.bounds, "transform": src.transform, "width": src.width,
                     "height": src.height, "nodata": src.nodata})
        except Exception as e:
            print(f"    警告: 读取元数据失败 for {tile_path}: {e}")
            continue

    if not tile_info_list:
        print("  没有可用的瓦片信息，镶嵌中止。")
        return False

    dst_bounds_left, dst_bounds_bottom, dst_bounds_right, dst_bounds_top = output_bounds
    output_transform = rasterio.Affine(first_tile_meta['transform'].a, 0.0, dst_bounds_left, 0.0,
                                       first_tile_meta['transform'].e, dst_bounds_top)
    output_width = int(round((dst_bounds_right - dst_bounds_left) / abs(output_transform.a)))
    output_height = int(round((dst_bounds_top - dst_bounds_bottom) / abs(output_transform.e)))

    if output_width <= 0 or output_height <= 0:
        print("  错误：计算得到的输出图像尺寸无效。")
        return False

    try:
        sum_pixels = np.zeros((output_height, output_width), dtype=np.float64)
        sum_weights = np.zeros((output_height, output_width), dtype=np.float32)
    except MemoryError:
        print("  错误: 内存分配失败，图像可能太大。")
        return False

    output_meta = first_tile_meta.copy()
    output_meta.update(
        {"driver": "GTiff", "height": output_height, "width": output_width, "transform": output_transform,
         "compress": "lzw"})
    nodata_val = first_tile_meta.get('nodata')

    for tile_info in tqdm(tile_info_list, desc="  镶嵌瓦片中"):
        with rasterio.open(tile_info['path']) as src_tile:
            tile_data = src_tile.read(1).astype(np.float64)
            tile_nodata = tile_info.get('nodata')
            is_nodata_mask = None
            if tile_nodata is not None: is_nodata_mask = (tile_data == tile_nodata); tile_data[is_nodata_mask] = np.nan
            weights = create_feather_weights(tile_info['height'], tile_info['width'], feather_percent)
            if is_nodata_mask is not None: weights[is_nodata_mask] = 0.0
            col_offset = int(round((tile_info['transform'].c - output_transform.c) / output_transform.a));
            row_offset = int(round((tile_info['transform'].f - output_transform.f) / output_transform.e))
            src_r_off, src_c_off = 0, 0;
            read_h, read_w = tile_info['height'], tile_info['width'];
            dst_r_off, dst_c_off = row_offset, col_offset;
            write_h, write_w = read_h, read_w
            if dst_c_off < 0: src_c_off = -dst_c_off; read_w += dst_c_off; write_w += dst_c_off; dst_c_off = 0
            if dst_r_off < 0: src_r_off = -dst_r_off; read_h += dst_r_off; write_h += dst_r_off; dst_r_off = 0
            if dst_c_off + write_w > output_width: overlap = (
                                                                         dst_c_off + write_w) - output_width; write_w -= overlap; read_w -= overlap
            if dst_r_off + write_h > output_height: overlap = (
                                                                          dst_r_off + write_h) - output_height; write_h -= overlap; read_h -= overlap
            if read_w <= 0 or read_h <= 0: continue
            valid_tile_data = tile_data[src_r_off:src_r_off + read_h, src_c_off:src_c_off + read_w];
            valid_weights = weights[src_r_off:src_r_off + read_h, src_c_off:src_c_off + read_w]
            target_sum_slice = sum_pixels[dst_r_off:dst_r_off + write_h, dst_c_off:dst_c_off + write_w];
            target_weight_slice = sum_weights[dst_r_off:dst_r_off + write_h, dst_c_off:dst_c_off + write_w]
            target_sum_slice += np.nan_to_num(valid_tile_data * valid_weights, nan=0.0);
            target_weight_slice += valid_weights

    with np.errstate(divide='ignore', invalid='ignore'):
        final_mosaic_data = sum_pixels / sum_weights
    if nodata_val is not None:
        final_mosaic_data[sum_weights == 0] = nodata_val
        final_mosaic_data = np.nan_to_num(final_mosaic_data, nan=nodata_val)
    else:
        final_mosaic_data[sum_weights == 0] = 0
        final_mosaic_data = np.nan_to_num(final_mosaic_data, nan=0)
    final_mosaic_data = final_mosaic_data.astype(first_tile_meta['dtype'])

    with rasterio.open(output_path, 'w', **output_meta) as dst:
        dst.write(final_mosaic_data, 1)
    print(f"  格网镶嵌完成！结果保存在: {output_path}")
    return True


if __name__ == "__main__":

    pixel_resolution = None
    if TILE_INDEX_PATH.exists():
        print(f"发现已存在的索引文件: {TILE_INDEX_PATH}")
        try:
            tiles_gdf = gpd.read_file(TILE_INDEX_PATH)
            print(f"索引加载成功，共 {len(tiles_gdf)} 个瓦片。")
            # 从加载的索引中获取像素分辨率
            first_tile_path = tiles_gdf['path'].iloc[0]
            with rasterio.open(first_tile_path) as src:
                pixel_resolution = abs(src.transform.a)
        except Exception as e:
            print(f"错误: 无法加载索引文件，将重新创建。错误信息: {e}")
            if TILE_INDEX_PATH.exists(): TILE_INDEX_PATH.unlink()

    if not TILE_INDEX_PATH.exists():
        print(f"未找到索引文件，正在从所有源目录创建新的索引...")
        all_tile_paths = []
        for input_dir in INPUT_DSM_TILES_DIRS:
            print(f"  正在扫描: {input_dir}")
            found_tiles = list(input_dir.glob("*.tif"))
            all_tile_paths.extend(found_tiles)
            print(f"  找到 {len(found_tiles)} 个文件。")

        if not all_tile_paths:
            print(f"错误: 在所有指定目录中均未找到.tif文件。")
            exit()

        print(f"总共找到 {len(all_tile_paths)} 个瓦片文件。")

        tile_index_data = []
        tile_crs = None
        for tile_path in tqdm(all_tile_paths, desc="构建索引"):
            try:
                with rasterio.open(tile_path) as src:
                    if tile_crs is None:
                        tile_crs = src.crs
                        pixel_resolution = abs(src.transform.a)  # 获取像元大小
                    geom = box(*src.bounds)
                    tile_index_data.append({'path': str(tile_path), 'geometry': geom})
            except Exception as e:
                print(f"警告: 无法读取 {tile_path}, 跳过. 错误: {e}")

        print("正在生成GeoDataFrame...")
        tiles_gdf = gpd.GeoDataFrame(tile_index_data, crs=tile_crs)

        print(f"索引创建完毕，共 {len(tiles_gdf)} 个瓦片。")
        print(f"正在将索引保存到: {TILE_INDEX_PATH} 以便下次使用...")
        try:
            tiles_gdf.to_file(TILE_INDEX_PATH, driver='GPKG')
            print("索引保存成功。")
        except Exception as e:
            print(f"警告: 保存索引文件失败！下次运行时将需要重新创建。错误: {e}")

    if pixel_resolution is None:
        print("错误：无法确定瓦片的像素分辨率。")
        exit()
      
    print(f"\n正在读取格网文件: {GRID_SHP_PATH}")
    try:
        grids_gdf = gpd.read_file(GRID_SHP_PATH)
        if grids_gdf.crs != tiles_gdf.crs:
            print(f"  坐标系不匹配，正在将格网重投影至瓦片坐标系...")
            grids_gdf = grids_gdf.to_crs(tiles_gdf.crs)
    except Exception as e:
        print(f"错误: 无法读取格网文件 {GRID_SHP_PATH}: {e}")
        exit()

    overlap_distance_map_units = OVERLAP_IN_PIXELS * pixel_resolution
    print(f"\n配置的像素重叠为 {OVERLAP_IN_PIXELS} px, 对应地图单位距离为: {overlap_distance_map_units:.2f}")
    print("-" * 50)

    # 使用 iterrows() 遍历所有格网
    for index, grid in tqdm(grids_gdf.iterrows(), total=len(grids_gdf), desc="处理所有格网"):
        grid_id = grid.get('Id', index)
        grid_geom = grid.geometry

        buffered_geom = grid_geom.buffer(overlap_distance_map_units)
        buffered_bounds = buffered_geom.bounds

        print(f"\n--- 开始处理格网 ID: {grid_id} ---")

        possible_matches_index = tiles_gdf.sindex.query(buffered_geom, predicate='intersects')
        possible_matches = tiles_gdf.iloc[possible_matches_index]
        intersecting_tiles_gdf = possible_matches[possible_matches.intersects(buffered_geom)]

        intersecting_files = intersecting_tiles_gdf['path'].tolist()

        if intersecting_files:
            output_mosaic_path = OUTPUT_MOSAIC_DIR / f"mosaic_grid_{grid_id}.tif"

            # 如果文件已存在，可以选择跳过
            if output_mosaic_path.exists():
                print(f"  文件 {output_mosaic_path.name} 已存在，跳过。")
                continue

            # 调用羽化镶嵌函数
            mosaic_tiles_for_grid(
                tile_files=intersecting_files,
                output_path=output_mosaic_path,
                output_bounds=buffered_bounds,  # 使用扩展后的边界
                feather_percent=FEATHER_WIDTH_PERCENT
            )
        else:
            print(f"  警告: 未找到任何与格网 {grid_id} 相交的瓦片。")

    print("\n" + "=" * 50)
    print(f"镶嵌结果保存在: {OUTPUT_MOSAIC_DIR}")
