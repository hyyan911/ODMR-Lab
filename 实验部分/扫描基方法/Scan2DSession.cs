﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ODMR_Lab.实验部分.扫描基方法.扫描范围;

namespace ODMR_Lab.实验部分.扫描基方法
{
    internal class Scan2DSession<T1, T2> : ScanHelper
    {
        /// <summary>
        /// 目标源1
        /// </summary>
        public T1 ScanSource1 { get; set; }

        /// <summary>
        /// 目标源2
        /// </summary>
        public T2 ScanSource2 { get; set; }

        /// <summary>
        /// 进度条设置操作（设备，当前进度）
        /// </summary>
        public Action<T1, T2, double> ProgressBarMethod = null;

        /// <summary>
        /// 设置状态操作（设备1,设备2,轴1位置,轴2位置）
        /// </summary>
        public Action<T1, T2, double, double> SetStateMethod = null;

        /// <summary>
        /// 扫描第一个点时进行的操作(源设备,轴1范围,轴2范围,轴1值,轴2值,输入参数,返回经过处理后的输入参数)
        /// </summary>
        public Func<T1, T2, D1NumricLinearScanRange, D1NumricLinearScanRange, double, double, List<object>, List<object>> FirstScanEvent = null;

        /// <summary>
        /// 正向扫描其他点时进行的操作(源设备,轴1范围,轴2范围,轴1值,轴2值,输入参数,返回经过处理后的输入参数)
        /// </summary>
        public Func<T1, T2, D1NumricLinearScanRange, D1NumricLinearScanRange, double, double, List<object>, List<object>> ScanEvent = null;
        /// <summary>
        /// 反向扫描其他点时进行的操作(源设备,轴1范围,轴2范围,轴1值,轴2值,输入参数,返回经过处理后的输入参数)
        /// </summary>
        public Func<T1, T2, D1NumricLinearScanRange, D1NumricLinearScanRange, double, double, List<object>, List<object>> ReverseScanEvent = null;

        /// <summary>
        /// 状态判断事件,此事件用来决定是否退出扫描步骤
        /// </summary>
        public Action StateJudgeEvent = null;


        /// <summary>
        /// 开始扫描(阻塞),ScanEvent输入参数为ps,返回值为最后一次执行ScanEvent得到的输出
        /// </summary>
        /// <param name="Scandir"></param>
        /// <param name="Lo"></param>
        /// <param name="Hi"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public List<object> BeginScan(D1NumricLinearScanRange axis1range, D1NumricLinearScanRange axis2range, double progressLo = 0, double progressHi = 100, params object[] ps)
        {
            SetProgress(progressLo);
            var scan1list = axis1range.ScanPoints;

            var scan2list = axis2range.ScanPoints;
            var scan2revlist = scan2list.ToArray().ToList();
            scan2revlist.Reverse();

            var progress = new D1NumricLinearScanRange(progressLo, progressHi, scan1list.Count * scan2list.Count * 2).ScanPoints;


            bool IsFirstScan = true;

            List<object> result = null;

            try
            {
                int ind = 0;
                bool IsReverse = false;
                for (int i = 0; i < scan1list.Count; i++)
                {
                    //正扫
                    for (int j = 0; j < scan2list.Count; j++)
                    {
                        //进行操作
                        if (IsFirstScan == true)
                        {
                            result = FirstScanEvent?.Invoke(ScanSource1, ScanSource2, axis1range, axis2range, scan1list[i], scan2list[j], ps.ToList());
                            StateJudgeEvent?.Invoke();
                            IsFirstScan = false;
                            ++ind;
                            continue;
                        }
                        result = ScanEvent?.Invoke(ScanSource1, ScanSource2, axis1range, axis2range, scan1list[i], scan2list[j], result);
                        StateJudgeEvent?.Invoke();
                        SetProgress(progress[ind]);
                        SetStateMethod?.Invoke(ScanSource1, ScanSource2, scan1list[i], scan2list[j]);
                        ++ind;
                    }
                    //反扫
                    for (int j = 0; j < scan2revlist.Count; j++)
                    {
                        //进行操作
                        if (IsFirstScan == true)
                        {
                            result = FirstScanEvent?.Invoke(ScanSource1, ScanSource2, axis1range, axis2range, scan1list[i], scan2list[j], ps.ToList());
                            StateJudgeEvent?.Invoke();
                            IsFirstScan = false;
                            ++ind;
                            continue;
                        }
                        result = ReverseScanEvent?.Invoke(ScanSource1, ScanSource2, axis1range, axis2range, scan1list[i], scan2list[j], result);
                        StateJudgeEvent?.Invoke();
                        SetProgress(progress[ind]);
                        SetStateMethod?.Invoke(ScanSource1, ScanSource2, scan1list[i], scan2list[j]);
                        ++ind;
                    }
                }
                return result;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        private void SetProgress(double v)
        {
            ProgressBarMethod?.Invoke(ScanSource1, ScanSource2, v);
        }
    }
}
