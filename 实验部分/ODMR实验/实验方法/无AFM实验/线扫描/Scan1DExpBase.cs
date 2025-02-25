﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ODMR_Lab.Python.LbviewHandler;
using System.Windows.Threading;
using ODMR_Lab.实验部分.磁场调节;
using System.Threading;
using System.Windows;
using ODMR_Lab.设备部分.位移台部分;
using ODMR_Lab.设备部分;
using ODMR_Lab.ODMR实验;
using HardWares.端口基类;
using ODMR_Lab.实验部分.扫描基方法;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.线扫描
{
    public abstract class Scan1DExpBase<T> : ODMRExperimentWithoutAFM
    {

        Scan1DSession<T> D1Session = new Scan1DSession<T>();

        public override void ODMRExpWithoutAFM()
        {
            T dev = GetScanSource();
            var range = GetScanRange();
            D1Session.FirstScanEvent = FirstScanEvent;
            D1Session.ScanEvent = ScanEvent;
            D1Session.ScanSource = dev;
            D1Session.ProgressBarMethod = new Action<T, double>((devi, val) =>
            {
                SetProgress(val);
            });
            D1Session.SetStateMethod = new Action<T, double>((devi, val1) =>
            {
                SetExpState(CreateThreadState(devi, val1));
            });
            D1Session.StateJudgeEvent = JudgeThreadEndOrResume;
            D1Session.BeginScan(range, 0, 100);
        }

        /// <summary>
        /// 设置扫描范围
        /// </summary>
        /// <returns></returns>
        public abstract ScanRange GetScanRange();

        /// <summary>
        /// 设置扫描源
        /// </summary>
        /// <returns></returns>
        public abstract T GetScanSource();

        /// <summary>
        /// 设置扫描范围
        /// </summary>
        /// <returns></returns>
        public abstract string CreateThreadState(T dev, double currentloc);

        /// <summary>
        /// 扫描操作，inputParams为空
        /// </summary>
        /// <param name="device"></param>
        /// <param name="locvalue"></param>
        /// <param name="inputParams"></param>
        /// <returns></returns>
        public abstract List<object> ScanEvent(T device, ScanRange range, double locvalue, List<object> inputParams);

        /// <summary>
        /// 第一次扫描操作，inputParams为空
        /// </summary>
        /// <param name="device"></param>
        /// <param name="locvalue"></param>
        /// <param name="inputParams"></param>
        /// <returns></returns>
        public abstract List<object> FirstScanEvent(T device, ScanRange range, double locvalue, List<object> inputParams);
    }
}
