﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ODMR_Lab.IO操作;

namespace ODMR_Lab.实验部分.磁场调节
{
    public class MagnetScanExpParams : ExpParamBase
    {
        #region 设备信息
        public Param<string> XDeviceName { get; set; } = new Param<string>("X方向设备", "", "XDeviceName");

        public Param<string> YDeviceName { get; set; } = new Param<string>("Y方向设备", "", "YDeviceName");

        public Param<string> ZDeviceName { get; set; } = new Param<string>("Z方向设备", "", "ZDeviceName");

        public Param<string> RotName { get; set; } = new Param<string>("旋转台设备", "", "RotName");
        #endregion

        #region 运行时参数

        /// <summary>
        /// X扫描结果
        /// </summary>
        public Param<double> XLoc { get; set; } = new Param<double>("X方向磁场峰值位置(mm)", double.NaN, "XLoc");

        /// <summary>
        /// Y扫描结果
        /// </summary>
        public Param<double> YLoc { get; set; } = new Param<double>("Y方向磁场峰值位置(mm)", double.NaN, "YLoc");

        /// <summary>
        /// Z扫描结果
        /// </summary>
        public Param<double> ZLoc { get; set; } = new Param<double>("Z方向位置(mm)", double.NaN, "ZLoc");

        /// <summary>
        /// Z扫描结果
        /// </summary>
        public Param<double> ZDistance { get; set; } = new Param<double>("Z方向位置和磁铁的距离(mm)", double.NaN, "ZDistance");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> Theta1 { get; set; } = new Param<double>("方位角θ1(度)", double.NaN, "Theta1");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> Theta2 { get; set; } = new Param<double>("方位角θ2(度)", double.NaN, "Theta2");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> Phi1 { get; set; } = new Param<double>("方位角φ1(度)", double.NaN, "Phi1");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> Phi2 { get; set; } = new Param<double>("方位角φ2(度)", double.NaN, "Phi2");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> CheckedPhi { get; set; } = new Param<double>("最终方位角φ(度)", double.NaN, "CheckedPhi");

        /// <summary>
        /// 角度扫描结果
        /// </summary>
        public Param<double> CheckedTheta { get; set; } = new Param<double>("最终方位角θ(度)", double.NaN, "CheckedTheta");
        #endregion
    }
}
