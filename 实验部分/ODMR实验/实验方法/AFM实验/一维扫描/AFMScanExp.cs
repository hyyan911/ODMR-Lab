﻿using Controls.Charts;
using Controls.Windows;
using HardWares.Lock_In;
using ODMR_Lab.IO操作;
using ODMR_Lab.ODMR实验;
using ODMR_Lab.基本控件;
using ODMR_Lab.基本控件.一维图表;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM;
using ODMR_Lab.实验部分.扫描基方法;
using ODMR_Lab.实验部分.扫描基方法.扫描范围;
using ODMR_Lab.设备部分;
using ODMR_Lab.设备部分.位移台部分;
using ODMR_Lab.设备部分.射频源_锁相放大器;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.AFM
{
    /// <summary>
    /// 一维扫描实验
    /// </summary>
    public class AFMScan1DExp : ODMRExperimentWithAFM
    {
        CustomScan1DLineSession<NanoStageInfo, NanoStageInfo> PointsScanSession = new CustomScan1DLineSession<NanoStageInfo, NanoStageInfo>();

        public override string ODMRExperimentName { get; set; } = "";
        public override string ODMRExperimentGroupName { get; set; } = "AFM线扫描";
        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<bool>("显示子实验窗口",true,"ShowSubMenu"),
            new Param<bool>("存储单点实验数据",true,"SaveSingleExpData"),
        };
        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>()
        {

        };

        public override List<KeyValuePair<DeviceTypes, Param<string>>> DeviceList { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>()
        {
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.AFM扫描台,new Param<string>("扫描台X","","ScannerX")),
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.AFM扫描台,new Param<string>("扫描台Y","","ScannerY")),
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.锁相放大器,new Param<string>("Lock In","","LockIn")),
        };
        public override List<ODMRExpObject> SubExperiments { get; set; } = new List<ODMRExpObject>();

        protected override List<KeyValuePair<string, Action>> AddInteractiveButtons()
        {
            return new List<KeyValuePair<string, Action>>();
        }
        public override List<ChartData1D> D1ChartDatas { get; set; } = new List<ChartData1D>();
        public override List<FittedData1D> D1FitDatas { get; set; } = new List<FittedData1D>();
        public override List<ChartData2D> D2ChartDatas { get; set; } = new List<ChartData2D>();

        public override bool IsAFMSubExperiment { get; protected set; } = false;

        public override void PreExpEventBeforeDropWithAFM()
        {
            LockinInfo info = GetDeviceByName("LockIn") as LockinInfo;
            //下针信息确认
            bool iscontinue = true;
            App.Current.Dispatcher.Invoke(() =>
            {
                if (MessageWindow.ShowMessageBox("下针信息确认", "当前振幅:" + info.Device.DemodR.ToString() + "\n" + "设定点:" + info.Device.SetPoint.ToString() + "\n"
                    + "P:" + info.Device.P.ToString() + "\n" + "I:" + info.Device.I.ToString() + "\n" + "D:" + info.Device.D.ToString() + "\n" + "是否继续?", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) != MessageBoxResult.Yes)
                {
                    iscontinue = false;
                }
            });
            if (!iscontinue) throw new Exception("实验已停止");
            //将位移台复位
            SetExpState("正在将位移台复位到零点...");
            NanoStageInfo infox = GetDeviceByName("ScannerX") as NanoStageInfo;
            NanoStageInfo infoy = GetDeviceByName("ScannerY") as NanoStageInfo;
            infox.Device.MoveToAndWait(D1ScanRange.StartPoint.X, 120000);
            infoy.Device.MoveToAndWait(D1ScanRange.StartPoint.Y, 120000);
        }

        public override void ODMRExpWithAFM()
        {
            NanoStageInfo dev1 = GetDeviceByName("ScannerX") as NanoStageInfo;
            NanoStageInfo dev2 = GetDeviceByName("ScannerY") as NanoStageInfo;

            PointsScanSession.FirstScanEvent = ScanEvent;
            PointsScanSession.ScanEvent = ScanEvent;
            PointsScanSession.ScanSource1 = dev1;
            PointsScanSession.ScanSource2 = dev2;
            PointsScanSession.ProgressBarMethod = new Action<NanoStageInfo, NanoStageInfo, double>((devx, devy, val) =>
             {
                 SetProgress(val);
             });
            PointsScanSession.SetStateMethod = new Action<NanoStageInfo, NanoStageInfo, Point>((devx, devy, v) =>
             {
                 SetExpState("当前扫描位置: X:" + Math.Round(v.X, 5).ToString() + " Y: " + Math.Round(v.Y, 5).ToString());
             });
            PointsScanSession.StateJudgeEvent = JudgeThreadEndOrResumeAction;
            PointsScanSession.BeginScan(D1ScanRange, 0, 100);
        }

        /// <summary>
        /// 扫描操作，inputParams为空
        /// </summary>
        /// <param name="device"></param>
        /// <param name="locvalue"></param>
        /// <param name="inputParams"></param>
        /// <returns></returns>
        public List<object> ScanEvent(NanoStageInfo scannerx, NanoStageInfo scannery, D1PointsScanRangeBase scanPoints, Point currentloc, List<object> inputParams)
        {
            SetExpState("正在移动到目标位置 X: " + Math.Round(currentloc.X, 5).ToString() + " Y: " + Math.Round(currentloc.Y, 5));
            scannerx.Device.MoveToAndWait(currentloc.X, 120000);
            scannery.Device.MoveToAndWait(currentloc.Y, 120000);

            SetExpState("正在进行实验...");
            //进行实验
            ODMRExpObject exp = RunSubExperimentBlock(0, GetInputParamValueByName("ShowSubMenu"));
            JudgeThreadEndOrResumeAction?.Invoke();
            //获取输出参数
            double value = 0;
            (Get1DChartDataSource("扫描点序号", "一维扫描数据")).Add(scanPoints.GetNearestIndex(currentloc));
            foreach (var item in exp.OutputParams)
            {
                if (item.RawValue is bool)
                {
                    value = (bool)item.RawValue ? 1 : 0;
                }
                if (item.RawValue is int)
                {
                    value = (int)item.RawValue;
                }
                if (item.RawValue is double)
                {
                    value = (double)item.RawValue;
                }
                if (item.RawValue is float)
                {
                    value = (float)item.RawValue;
                }
                var data = Get1DChartData(ODMRExperimentGroupName + ":" + ODMRExperimentName + ":" + item.Description, "一维扫描数据");
                if (data == null)
                {
                    data = new NumricChartData1D(ODMRExperimentGroupName + ":" + ODMRExperimentName + ":" + item.Description, "一维扫描数据", ChartDataType.Y);
                    D1ChartDatas.Add(data);
                    UpdatePlotChart();
                }
                (data as NumricChartData1D).Data.Add(value);
            }
            //刷新形貌数据
            var afm = Get1DChartDataSource("AFM形貌数据(PID输出)", "一维扫描数据");
            afm.Add((GetDeviceByName("LockIn") as LockinInfo).Device.PIDValue);
            UpdatePlotChartFlow(true);

            //设置一维图表
            if (exp is ODMRExperimentWithoutAFM)
            {
                var plotpacks = (exp as ODMRExperimentWithoutAFM).GetD1PlotPacks();
                foreach (var plot in plotpacks)
                {
                    if (!plot.IsContinusPlot)
                    {
                        if (Get1DChartData(plot.Name, plot.GroupName) == null)
                        {
                            D1ChartDatas.Add(new NumricChartData1D(plot.Name, plot.GroupName, plot.AxisType) { Data = plot.ChartData });
                        }
                    }
                    else
                    {
                        D1ChartDatas.Add(new NumricChartData1D("X:" + Math.Round(currentloc.X, 5).ToString() + "Y:" + Math.Round(currentloc.Y, 5).ToString() + " " + plot.Name, plot.GroupName, plot.AxisType) { Data = plot.ChartData });
                    }
                }
            }

            UpdatePlotChartFlow(true);
            return new List<object>();
        }

        public override void AfterExpEventWithAFM()
        {

        }

        public override void PreExpEventWithAFM()
        {
            D1ChartDatas.Clear();
            D1ChartDatas.Add(new NumricChartData1D("扫描点序号", "一维扫描数据", ChartDataType.X));
            D1ChartDatas.Add(new NumricChartData1D("AFM形貌数据", "一维扫描数据", ChartDataType.Y));
            UpdatePlotChart();
        }

        public override LockinInfo GetLockIn()
        {
            return GetDeviceByName("LockIn") as LockinInfo;
        }

        public override bool PreConfirmProcedure()
        {
            if (D1ScanRange == null)
            {
                MessageWindow.ShowTipWindow("扫描范围未设置", Window.GetWindow(ParentPage));
                return false;
            }

            if (MessageWindow.ShowMessageBox("提示", "是否要继续?此操作将清除原先的实验数据,并且将执行下针操作", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) != MessageBoxResult.Yes)
            {
                return false;
            }

            //扫描范围信息
            if (MessageWindow.ShowMessageBox("扫描范围确认", "扫描范围类型:" + D1ScanRange.ScanName + "\n" + D1ScanRange.GetDescription() + "\n" + "是否继续?", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) != MessageBoxResult.Yes)
            {
                return false;
            }

            return true;
        }
    }
}
