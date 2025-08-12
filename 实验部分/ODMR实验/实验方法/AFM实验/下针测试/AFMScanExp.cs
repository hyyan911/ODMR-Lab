﻿using Controls.Charts;
using Controls.Windows;
using HardWares.Lock_In;
using ODMR_Lab.IO操作;
using ODMR_Lab.ODMR实验;
using ODMR_Lab.基本控件;
using ODMR_Lab.基本控件.一维图表;
using ODMR_Lab.基本窗口;
using ODMR_Lab.实验部分.ODMR实验.实验方法.ScanCore;
using ODMR_Lab.实验部分.ODMR实验.实验方法.其他;
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
    public class AFMScanTestExp : ODMRExperimentWithAFM
    {
        public override bool Is1DScanExp { get; set; } = false;
        public override bool Is2DScanExp { get; set; } = false;

        public override string ODMRExperimentName { get; set; } = "AFM下针测试";
        public override string ODMRExperimentGroupName { get; set; } = "AFM下针测试";
        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<bool>("尝试多次下针",true,"MultiAFMDrop"),
            new Param<int>("最大尝试下针次数",5,"DropCount"),
            new Param<double>("尝试下针失败后样品移动量",0.005,"DropDet"),
            new Param<bool>("样品轴升高方向反向",false,"SampleAxisReverse"),
            new Param<bool>("是否悬浮测量",false,"IsFloatScan"),
            new Param<double>("悬浮测量距离(V)",0.1,"FloatHeight"),
            new Param<double>("最大限制电压(V)",10,"UpperLimit"),
        };
        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>()
        {

        };

        public override List<KeyValuePair<DeviceTypes, Param<string>>> DeviceList { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>()
        {
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.锁相放大器,new Param<string>("Lock In","","LockIn")),
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.样品位移台,new Param<string>("样品Z轴","","SampleZ")),
        };
        public override List<ODMRExpObject> SubExperiments { get; set; } = new List<ODMRExpObject>()
        {
        };

        protected override List<KeyValuePair<string, Action>> AddInteractiveButtons()
        {
            List<KeyValuePair<string, Action>> btns = new List<KeyValuePair<string, Action>>();
            btns.Add(new KeyValuePair<string, Action>("移动扫描台到指定位置", MoveScanner));
            return btns;
        }
        public override List<ChartData1D> D1ChartDatas { get; set; } = new List<ChartData1D>();
        public override List<FittedData1D> D1FitDatas { get; set; } = new List<FittedData1D>();
        public override List<ChartData2D> D2ChartDatas { get; set; } = new List<ChartData2D>();

        public override bool IsAFMSubExperiment { get; protected set; } = false;

        public override void PreExpEventBeforeDropWithAFM()
        {
        }

        public override void ODMRExpWithAFM()
        {
            //悬浮测量
            if (GetInputParamValueByName("IsFloatScan"))
            {
                AFMFloatDrop drop = new AFMFloatDrop();
                var result = drop.CoreMethod(new List<object>() { GetInputParamValueByName("UpperLimit"), GetInputParamValueByName("FloatHeight") }, GetDeviceByName("LockIn"));
                if ((bool)result[0] == false)
                {
                    throw new Exception("下针失败，实验已结束");
                }
            }
            //等待实验结束
            SetExpState("下针完成,等待手动结束并撤针...");
            while (true)
            {
                JudgeThreadEndOrResumeAction?.Invoke();
            }
        }

        public override void AfterExpEventWithAFM()
        {
            if (GetInputParamValueByName("IsFloatScan"))
            {
                (GetDeviceByName("LockIn") as LockinInfo).Device.PIDOutputUpperLimit = GetInputParamValueByName("UpperLimit");
            }
        }

        public override void PreExpEventWithAFM()
        {
        }

        public override LockinInfo GetLockIn()
        {
            return GetDeviceByName("LockIn") as LockinInfo;
        }

        public override bool PreConfirmProcedure()
        {
            if (MessageWindow.ShowMessageBox("提示", "是否要继续?此操作将清除原先的实验数据,并且将执行下针操作", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) != MessageBoxResult.Yes)
            {
                return false;
            }
            return true;
        }

        protected override bool GetAllowMultiDrop()
        {
            return GetInputParamValueByName("MultiAFMDrop");
        }

        protected override int GetMaxDropCount()
        {
            return GetInputParamValueByName("DropCount");
        }

        protected override double GetDropDet()
        {
            return GetInputParamValueByName("DropDet");
        }

        protected override bool GetSampleAxisReverse()
        {
            return GetInputParamValueByName("SampleAxisReverse");
        }
        protected override NanoStageInfo GetSampleZ()
        {
            return GetDeviceByName("SampleZ") as NanoStageInfo;
        }
    }
}
