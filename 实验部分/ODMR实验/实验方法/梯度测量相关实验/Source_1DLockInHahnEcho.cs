﻿using Controls.Windows;
using ODMR_Lab.IO操作;
using ODMR_Lab.ODMR实验;
using ODMR_Lab.基本控件;
using ODMR_Lab.基本控件.一维图表;
using ODMR_Lab.实验部分.ODMR实验.实验方法.AFM;
using ODMR_Lab.实验部分.ODMR实验.实验方法.其他;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.点实验.脉冲实验;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM实验;
using ODMR_Lab.实验部分.扫描基方法.扫描范围;
using ODMR_Lab.实验部分.扫描基方法;
using ODMR_Lab.设备部分;
using ODMR_Lab.设备部分.射频源_锁相放大器;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using ODMR_Lab.设备部分.源表;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.梯度测量相关实验
{
    internal class Source_1DLockInHahnEcho : VibrationExpBase
    {
        public override List<KeyValuePair<DeviceTypes, Param<string>>> PulseExpDevices { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>()
        {
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.锁相放大器,new Param<string>("Lock In","","LockIn")),
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.源表,new Param<string>("电压源","","VoltSource")),
        };
        public override bool Is1DScanExp { get; set; } = true;
        public override bool Is2DScanExp { get; set; } = false;
        public override string ODMRExperimentName { get; set; } = "变样品电压一维扫描";
        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<double>("源表电压起始值（V）",0,"SignalStart"),
            new Param<double>("源表电压终止值（V）",1e-3,"SignalEnd"),
            new Param<int>("变压扫描点数",10,"SignalCount"),
            new Param<double>("源表电压变压步长（V）",1e-4,"SignalRampGap"),
        };
        public override List<ParamB> OutputParams { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public override List<ChartData1D> D1ChartDatas { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public override List<FittedData1D> D1FitDatas { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public override List<ChartData2D> D2ChartDatas { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public override void AfterExpEventWithoutAFM()
        {

        }

        public override List<ParentPlotDataPack> GetD1PlotPacks()
        {
            return new List<ParentPlotDataPack>();
        }

        public override void ODMRExpWithoutAFM()
        {
            D1NumricLinearScanRange range = new D1NumricLinearScanRange(GetInputParamValueByName("SignalStart"), GetInputParamValueByName("SignalEnd"), GetInputParamValueByName("SignalCount"));
            Scan1DSession<object> session = new Scan1DSession<object>();
            session.SetStateMethod = new Action<object, double>((obj, val) =>
            {
                SetExpState("当前电压(V):" + val.ToString());
            });
            session.ScanSource = new object();
            session.ProgressBarMethod = new Action<object, double>((obj, val) =>
            {
                SetProgress(val);
            });
            session.StateJudgeEvent = JudgeThreadEndOrResumeAction;
            session.FirstScanEvent = ScanEvent;
            session.ScanEvent = ScanEvent;
            session.BeginScan(range, 0, 100);
        }

        private List<object> ScanEvent(object arg1, D1NumricScanRangeBase @base, double arg3, List<object> list)
        {
            (GetDeviceByName("VoltSource") as PowerMeterInfo).Device.VoltageRampStep = GetInputParamValueByName("SignalRampGap");
            //设置驱动电压
            (GetDeviceByName("VoltSource") as PowerMeterInfo).Device.TargetVoltage = arg3;

            var exp = RunSubExperimentBlock(0, true);
            foreach (var item in exp.D1ChartDatas)
            {
                item.Name = "电压" + arg3.ToString() + "V --" + item.Name;
                D1ChartDatas.Add(item);
            }
            ;
            UpdatePlotChart();
            UpdatePlotChartFlow(true);
            JudgeThreadEndOrResumeAction?.Invoke();
            return new List<object>();
        }

        public override bool PreConfirmProcedure()
        {
            if (DropConfirm(GetDeviceByName("LockIn") as LockinInfo) == false) return false;
            return true;
        }

        public override void PreExpEventWithoutAFM()
        {
            //配置子实验
            (SubExperiments[0] as LockIn1DScan).D1ScanRange = D1ScanRange;
        }

        protected override List<KeyValuePair<string, Action>> AddPulseInteractiveButtons()
        {
            return new List<KeyValuePair<string, Action>>();
        }

        protected override List<ODMRExpObject> GetSubExperiments()
        {
            return new List<ODMRExpObject>()
            {
               new LockIn1DScan()
            };
        }
    }
}
